"""Tests for the fragment-spectrum endpoint, including a brute-force cross-check."""

from __future__ import annotations

import itertools
import math
import random
from fractions import Fraction

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def parse_fraction(text: str) -> Fraction:
    n, d = text.split("/")
    return Fraction(int(n), int(d))


def fmt(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def brute_force_fragment(payload: dict, fragment_counts, target_parent: int):
    """Enumerate every isotope assignment and every subset-selection path.

    Given the isotope values, each group contributes exactly
    ``fragment_counts[g]`` atoms drawn uniformly without replacement, so a
    subset of size k in a group of size n carries probability 1/C(n, k).
    Returns (parent-peak probability, {fragment shift: conditional prob}).
    """
    groups = payload["atom_groups"]
    flat = []
    group_atoms = []
    for group in groups:
        atoms = [
            (iso["shift_mda"], parse_fraction(iso["abundance"]))
            for iso in group["isotopes"]
        ]
        indices = list(range(len(flat), len(flat) + group["count"]))
        for _ in range(group["count"]):
            flat.append(atoms)
        group_atoms.append(indices)

    subset_prob = Fraction(
        1,
        math.prod(
            math.comb(len(indices), k)
            for indices, k in zip(group_atoms, fragment_counts)
        ),
    )

    parent_probability = Fraction(0)
    joint = {}
    for assignment in itertools.product(*flat):
        parent_shift = sum(shift for shift, _ in assignment)
        probability = Fraction(1)
        for _, abundance in assignment:
            probability *= abundance
        if parent_shift != target_parent:
            continue
        parent_probability += probability
        selections = [
            list(itertools.combinations(indices, k))
            for indices, k in zip(group_atoms, fragment_counts)
        ]
        for chosen in itertools.product(*selections):
            drawn = {index for subset in chosen for index in subset}
            fragment_shift = sum(assignment[i][0] for i in drawn)
            joint[fragment_shift] = (
                joint.get(fragment_shift, Fraction(0))
                + probability * subset_prob
            )

    if parent_probability == 0:
        return Fraction(0), {}
    conditional = {
        shift: weight / parent_probability
        for shift, weight in joint.items()
        if weight != 0
    }
    return parent_probability, conditional


def basic_payload():
    return {
        "atom_groups": [
            {
                "count": 3,
                "isotopes": [
                    {"shift_mda": 0, "abundance": "1/2"},
                    {"shift_mda": 1, "abundance": "1/2"},
                ],
            }
        ],
        "fragment_counts": [1],
        "target_parent_shift_mda": 1,
    }


def test_basic_single_group():
    # Parent shift 1 among 3 fair binary atoms: C(3,1)/8 = 3/8. Given exactly
    # one heavy atom, the one-in-three fragment draw finds it with prob 1/3,
    # so fragment shifts are 0 -> 2/3 and 1 -> 1/3.
    resp = client.post("/v1/fragment-spectrum", json=basic_payload())
    assert resp.status_code == 200
    assert resp.json() == {
        "reachable": True,
        "probability": "3/8",
        "peaks": [
            {"shift_mda": 0, "probability": "2/3"},
            {"shift_mda": 1, "probability": "1/3"},
        ],
    }


def test_draw_nothing_fragment_is_residual_only():
    payload = basic_payload()
    payload["fragment_counts"] = [0]
    payload["target_parent_shift_mda"] = 2
    resp = client.post("/v1/fragment-spectrum", json=payload)
    assert resp.status_code == 200
    assert resp.json() == {
        "reachable": True,
        "probability": "3/8",
        "peaks": [{"shift_mda": 0, "probability": "1/1"}],
    }


def test_draw_everything_fragment_is_the_parent():
    payload = basic_payload()
    payload["fragment_counts"] = [3]
    payload["target_parent_shift_mda"] = 2
    resp = client.post("/v1/fragment-spectrum", json=payload)
    assert resp.status_code == 200
    assert resp.json() == {
        "reachable": True,
        "probability": "3/8",
        "peaks": [{"shift_mda": 2, "probability": "1/1"}],
    }


def test_unreachable_parent_shift():
    payload = basic_payload()
    payload["target_parent_shift_mda"] = 9
    resp = client.post("/v1/fragment-spectrum", json=payload)
    assert resp.status_code == 200
    assert resp.json() == {"reachable": False, "probability": "0/1"}
    assert "peaks" not in resp.json()


def test_zero_abundance_isotope_never_drawn():
    payload = {
        "atom_groups": [
            {
                "count": 2,
                "isotopes": [
                    {"shift_mda": 0, "abundance": "1/1"},
                    {"shift_mda": 5, "abundance": "0/1"},
                ],
            }
        ],
        "fragment_counts": [1],
        "target_parent_shift_mda": 0,
    }
    resp = client.post("/v1/fragment-spectrum", json=payload)
    assert resp.status_code == 200
    assert resp.json() == {
        "reachable": True,
        "probability": "1/1",
        "peaks": [{"shift_mda": 0, "probability": "1/1"}],
    }


def test_multi_group_matches_brute_force():
    payload = {
        "atom_groups": [
            {
                "count": 3,
                "isotopes": [
                    {"shift_mda": 0, "abundance": "3/4"},
                    {"shift_mda": 1, "abundance": "1/4"},
                ],
            },
            {
                "count": 2,
                "isotopes": [
                    {"shift_mda": 0, "abundance": "2/3"},
                    {"shift_mda": 2, "abundance": "1/3"},
                ],
            },
        ]
    }
    max_shift = 3 * 1 + 2 * 2
    rng = random.Random(0)
    for _ in range(40):
        counts = [rng.randint(0, group["count"]) for group in payload["atom_groups"]]
        target = rng.randint(0, max_shift + 1)
        resp = client.post(
            "/v1/fragment-spectrum",
            json={
                **payload,
                "fragment_counts": counts,
                "target_parent_shift_mda": target,
            },
        )
        assert resp.status_code == 200
        body = resp.json()
        probability, conditional = brute_force_fragment(payload, counts, target)

        if probability == 0:
            assert body == {"reachable": False, "probability": "0/1"}
            continue

        assert body["reachable"] is True
        assert body["probability"] == fmt(probability)
        assert parse_fraction(body["probability"]) == probability

        peaks = body["peaks"]
        shifts = [peak["shift_mda"] for peak in peaks]
        # Non-zero peaks only, ascending by shift.
        assert shifts == sorted(shifts)
        assert all(parse_fraction(peak["probability"]) != 0 for peak in peaks)
        got = {
            peak["shift_mda"]: parse_fraction(peak["probability"]) for peak in peaks
        }
        assert got == conditional
        # The conditional distribution is normalized and reduced.
        assert sum(got.values()) == 1
        for peak in peaks:
            n, d = peak["probability"].split("/")
            assert Fraction(int(n), int(d)).denominator == int(d)
        # Fragment shifts never exceed the observed parent shift.
        assert max(shifts, default=0) <= target


def test_ternary_three_group_matches_brute_force():
    payload = {
        "atom_groups": [
            {
                "count": 2,
                "isotopes": [
                    {"shift_mda": 0, "abundance": "1/1"},
                    {"shift_mda": 5, "abundance": "0/1"},
                ],
            },
            {
                "count": 4,
                "isotopes": [
                    {"shift_mda": 0, "abundance": "7/10"},
                    {"shift_mda": 1, "abundance": "1/5"},
                    {"shift_mda": 3, "abundance": "1/10"},
                ],
            },
            {
                "count": 1,
                "isotopes": [
                    {"shift_mda": 0, "abundance": "1/2"},
                    {"shift_mda": 7, "abundance": "1/2"},
                ],
            },
        ]
    }
    max_shift = 2 * 5 + 4 * 3 + 7
    rng = random.Random(1)
    for _ in range(40):
        counts = [rng.randint(0, group["count"]) for group in payload["atom_groups"]]
        target = rng.randint(0, max_shift + 1)
        body = client.post(
            "/v1/fragment-spectrum",
            json={
                **payload,
                "fragment_counts": counts,
                "target_parent_shift_mda": target,
            },
        ).json()
        probability, conditional = brute_force_fragment(payload, counts, target)
        if probability == 0:
            assert body == {"reachable": False, "probability": "0/1"}
            continue
        assert body["probability"] == fmt(probability)
        got = {
            peak["shift_mda"]: parse_fraction(peak["probability"])
            for peak in body["peaks"]
        }
        assert got == conditional
        assert sum(got.values()) == 1


def test_parent_probability_matches_envelope_peak():
    payload = {
        "atom_groups": [
            {
                "count": 4,
                "isotopes": [
                    {"shift_mda": 0, "abundance": "3/5"},
                    {"shift_mda": 2, "abundance": "2/5"},
                ],
            },
            {
                "count": 2,
                "isotopes": [
                    {"shift_mda": 0, "abundance": "5/8"},
                    {"shift_mda": 1, "abundance": "3/8"},
                ],
            },
        ]
    }
    envelope = client.post("/v1/envelope", json={**payload, "k": 20}).json()
    probs = {peak["shift_mda"]: peak["probability"] for peak in envelope["peaks"]}
    for target in range(0, 4 * 2 + 2 * 1 + 1):
        body = client.post(
            "/v1/fragment-spectrum",
            json={
                **payload,
                "fragment_counts": [2, 1],
                "target_parent_shift_mda": target,
            },
        ).json()
        if target in probs:
            assert body["reachable"] is True
            assert body["probability"] == probs[target]
        else:
            assert body == {"reachable": False, "probability": "0/1"}


@pytest.mark.parametrize(
    "mutation",
    [
        # fragment_counts range / shape / type
        lambda p: p.update(fragment_counts=[4]),
        lambda p: p.update(fragment_counts=[-1]),
        lambda p: p.update(fragment_counts=[]),
        lambda p: p.update(fragment_counts=[1, 0]),
        lambda p: p.update(fragment_counts=["1"]),
        lambda p: p.update(fragment_counts=[1.0]),
        lambda p: p.update(fragment_counts=[True]),
        lambda p: p.pop("fragment_counts"),
        # target parent shift out of range / wrong type
        lambda p: p.update(target_parent_shift_mda=-1),
        lambda p: p.update(target_parent_shift_mda=2001),
        lambda p: p.update(target_parent_shift_mda=1.5),
        lambda p: p.update(target_parent_shift_mda="1"),
        lambda p: p.pop("target_parent_shift_mda"),
        # unknown field
        lambda p: p.update(extra=1),
        # atom_groups violations reused verbatim from /v1/envelope
        lambda p: p["atom_groups"][0]["isotopes"][1].update(abundance="2/3"),
        lambda p: p["atom_groups"][0]["isotopes"][1].update(shift_mda=0),
        lambda p: p["atom_groups"][0].update(count=0),
        lambda p: p.update(atom_groups=[]),
    ],
)
def test_invalid_payloads_return_422(mutation):
    payload = basic_payload()
    mutation(payload)
    resp = client.post("/v1/fragment-spectrum", json=payload)
    assert resp.status_code == 422


def test_malformed_json_returns_422():
    resp = client.post(
        "/v1/fragment-spectrum",
        content="{not json",
        headers={"Content-Type": "application/json"},
    )
    assert resp.status_code == 422
