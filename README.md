# 同位素峰型计算 API（Isotope Envelope API）

基于 FastAPI 的本地服务：给定若干原子组（每组有 2~3 种同位素及其质量偏移与丰度），
在**精确有理数**语义下计算完整的同位素峰型（envelope）。各原子独立地取同位素，
所有路径按质量偏移合并，全程使用 `fractions.Fraction` 运算，**不中途裁剪任何路径**。
此外还提供：给定单个观测峰的**归因**接口（该峰的精确概率与条件期望原子数），针对
观测峰概率区间求解多组分混合比例的接口，在固定采样预算下规划**采集窗口**的接口，
有限离子采样下多个（可重叠）检测窗口同时达标的**联合检出概率**接口，在共同
整数校准偏移下对 2～4 个互斥候选做精确贝叶斯**候选/偏移后验**推断的接口，以及在
两个等先验候选之间规划**判别窗口**、精确最大化判别准确率的接口，
同样全程精确有理数运算。

## 启动

一条命令（Docker Compose）：

```bash
docker compose up --build
```

服务监听 **8000** 端口。交互式文档见 `http://localhost:8000/docs`。

也可以不用 Docker 本地运行：

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## API

### `POST /v1/envelope`

请求体：

```json
{
  "atom_groups": [
    {
      "count": 2,
      "isotopes": [
        {"shift_mda": 0, "abundance": "3/4"},
        {"shift_mda": 1, "abundance": "1/4"}
      ]
    }
  ],
  "k": 2
}
```

字段与约束（任一不满足即返回 **422**）：

| 字段 | 约束 |
| --- | --- |
| `atom_groups` | 1～4 组 |
| `k` | 整数，1～20 |
| `count` | 整数，1～40；各组 `count` 总和 ≤ 80 |
| `isotopes` | 每组 2～3 个 |
| `shift_mda` | 非负整数；**组内不得重复** |
| `abundance` | 字符串 `"n/d"`，`0 ≤ n ≤ d`，`1 ≤ d ≤ 1000`；**每组丰度之和必须精确等于 1** |
| 整体 | 最大可能总偏移 `Σ count × max(shift_mda)` ≤ 2000 |

响应（200）：

```json
{
  "peaks": [
    {"shift_mda": 0, "probability": "9/16"},
    {"shift_mda": 1, "probability": "3/8"},
    {"shift_mda": 2, "probability": "1/16"}
  ],
  "top_peaks": [
    {"shift_mda": 0, "probability": "9/16"},
    {"shift_mda": 1, "probability": "3/8"}
  ]
}
```

- `peaks`：全部非零峰，按 `shift_mda` 升序；概率为**已约分**的 `"n/d"`。
- `top_peaks`：前 `k` 个峰，按概率降序，同概率按 `shift_mda` 升序。

### 调用示例

```bash
curl -s -X POST http://localhost:8000/v1/envelope \
  -H 'Content-Type: application/json' \
  -d '{
    "atom_groups": [
      {"count": 2, "isotopes": [
        {"shift_mda": 0, "abundance": "3/4"},
        {"shift_mda": 1, "abundance": "1/4"}
      ]}
    ],
    "k": 2
  }'
```

### `POST /v1/peak-attribution`

给定目标质量偏移，返回该偏移处的**峰概率**，以及在**已观测到此偏移**的条件下，
每组中各同位素的**原子个数期望**。所有能产生该偏移的取值路径都按其原概率参与条件化
（不重新归一化之外的任何裁剪），全程 `Fraction` 运算；零丰度同位素的期望为 `"0/1"`。

请求体：

```json
{
  "atom_groups": [
    {
      "count": 2,
      "isotopes": [
        {"shift_mda": 0, "abundance": "3/4"},
        {"shift_mda": 1, "abundance": "1/4"}
      ]
    }
  ],
  "target_shift_mda": 1
}
```

字段与约束（任一不满足即返回 **422**）：

| 字段 | 约束 |
| --- | --- |
| `atom_groups` | 与 `/v1/envelope` 完全相同的校验（1～4 组、各组 `count` 总和 ≤ 80、最大可能总偏移 ≤ 2000 等） |
| `target_shift_mda` | 整数，0～2000 |

响应（200）。目标峰可达时：

```json
{
  "reachable": true,
  "probability": "3/8",
  "expected_counts": [["1/1", "1/1"]]
}
```

目标峰概率为零时（没有任何路径能产生该偏移）：

```json
{"reachable": false, "probability": "0/1"}
```

- `probability`：目标偏移的（无条件）峰概率，已约分的 `"n/d"`。
- `expected_counts`：二维 `n/d` 数组，按**组顺序**与组内**同位素输入顺序**排列；
  `expected_counts[g][i]` 即在目标偏移条件下第 `g` 组第 `i` 种同位素的期望原子数；
  每组期望值之和等于该组的 `count`。仅在 `reachable=true` 时返回。
- 目标峰不可达时返回 `reachable=false`、`probability="0/1"`，不返回 `expected_counts`。

### 调用示例

```bash
curl -s -X POST http://localhost:8000/v1/peak-attribution \
  -H 'Content-Type: application/json' \
  -d '{
    "atom_groups": [
      {"count": 2, "isotopes": [
        {"shift_mda": 0, "abundance": "3/4"},
        {"shift_mda": 1, "abundance": "1/4"}
      ]}
    ],
    "target_shift_mda": 1
  }'
```

### `POST /v1/mixture-bounds`给定 2～4 个组分（各自复用完整峰型）和 1～8 条观测（某偏移处混合峰概率的闭区间），
在**精确有理数**语义下求解非负且总和为 1 的组分比例 `w`，使每条已测偏移处的
混合概率 `Σ w_i · P_i(shift)`（缺峰按 0）落在对应闭区间内。未测偏移不加约束，
也不重新归一化。

请求体：

```json
{
  "components": [
    {
      "atom_groups": [
        {"count": 1, "isotopes": [
          {"shift_mda": 0, "abundance": "1/2"},
          {"shift_mda": 1, "abundance": "1/2"}
        ]}
      ]
    },
    {
      "atom_groups": [
        {"count": 1, "isotopes": [
          {"shift_mda": 0, "abundance": "1/2"},
          {"shift_mda": 2, "abundance": "1/2"}
        ]}
      ]
    }
  ],
  "observations": [
    {"shift_mda": 1, "lower": "1/4", "upper": "1/2"}
  ]
}
```

字段与约束（任一不满足即返回 **422**）：

| 字段 | 约束 |
| --- | --- |
| `components` | 2～4 个组分；峰型相同的组分仍保留为独立组分 |
| `components[].atom_groups` | 与 `/v1/envelope` 完全相同的校验（1～4 组、各组 `count` 总和 ≤ 80、最大可能总偏移 ≤ 2000 等） |
| `observations` | 1～8 条 |
| `observations[].shift_mda` | 整数，0～2000；各观测之间**不得重复** |
| `observations[].lower` / `upper` | 与 `abundance` 相同的 `"n/d"` 格式及范围，且 `lower ≤ upper` |

响应（200）。无解时：

```json
{"feasible": false}
```

有解时：

```json
{
  "feasible": true,
  "ranges": [
    {"min": "1/2", "max": "1/1"},
    {"min": "0/1", "max": "1/2"}
  ],
  "weights": ["1/2", "1/2"]
}
```

- `ranges`：每个组分比例在**全部可行解**中的最小/最大值，按组分输入顺序排列。
- `weights`：按组分输入顺序字典序最小的可行比例向量。
- 所有分数均为已约分的 `"n/d"`。

### 调用示例

```bash
curl -s -X POST http://localhost:8000/v1/mixture-bounds \
  -H 'Content-Type: application/json' \
  -d '{
    "components": [
      {"atom_groups": [
        {"count": 1, "isotopes": [
          {"shift_mda": 0, "abundance": "1/2"},
          {"shift_mda": 1, "abundance": "1/2"}
        ]}
      ]},
      {"atom_groups": [
        {"count": 1, "isotopes": [
          {"shift_mda": 0, "abundance": "1/2"},
          {"shift_mda": 2, "abundance": "1/2"}
        ]}
      ]}
    ],
    "observations": [
      {"shift_mda": 1, "lower": "1/4", "upper": "1/2"}
    ]
  }'
```

### `GET /healthz`

健康检查，返回 `{"status": "ok"}`。

### `POST /v1/acquisition-plan`

在复用完整峰型（**不**先截取 `top_peaks`）的基础上，于 **0～2000 mDa 质量偏移轴**上规划
采集窗口：选择至多 `max_windows` 个**互不相交**的整数闭区间 `[l, r]`，使窗口并集所覆盖
的**峰概率之和**最大，全程精确有理数运算。

**窗口计费**：每个窗口 `[l, r]` 占用 `r - l + 1` 个采样格；窗口内没有峰的格子同样计费，
所有窗口的格数之和即 `total_bins`，必须满足 `total_bins ≤ budget`。例如峰在 0 和 10，
单窗口 `[0, 10]` 计费 11 格；改成 `[0, 0]` 与 `[10, 10]` 两个窗口则只计费 2 格。

并列时依次按以下顺序决胜：

1. 覆盖概率相同 → `total_bins` 更少者胜（等价地，最优方案的窗口端点必落在某个峰上）；
2. 仍并列 → 窗口数更少者胜；
3. 仍并列 → 窗口按起点升序排序、端点展平为 `[l1, r1, l2, r2, …]` 后，字典序更小者胜。

请求体：

```json
{
  "atom_groups": [
    {
      "count": 2,
      "isotopes": [
        {"shift_mda": 0, "abundance": "3/4"},
        {"shift_mda": 1, "abundance": "1/4"}
      ]
    }
  ],
  "max_windows": 2,
  "budget": 5
}
```

字段与约束（任一不满足即返回 **422**）：

| 字段 | 约束 |
| --- | --- |
| `atom_groups` | 与 `/v1/envelope` 完全相同的校验（1～4 组、各组 `count` 总和 ≤ 80、最大可能总偏移 ≤ 2000 等） |
| `max_windows` | 整数，1～6；窗口数**至多**为该值 |
| `budget` | 整数，1～200；`total_bins = Σ(r−l+1)` 不得超过该值 |

响应（200）：

```json
{
  "windows": [[0, 2]],
  "total_bins": 3,
  "captured_probability": "1/1"
}
```

- `windows`：`[l, r]` 数组，按起点 `l` 升序；两窗口互不相交。
- `total_bins`：实际计费格数 `Σ(r−l+1)`，保证 `≤ budget`。
- `captured_probability`：窗口并集内全部峰概率之和，已约分的 `"n/d"`。

上例中 `[0, 2]`（1 窗、3 格）与 `[0, 1] + [2, 2]`（2 窗、3 格）覆盖概率相同，决胜规则
选窗口更少的 `[[0, 2]]`。

### 调用示例

```bash
curl -s -X POST http://localhost:8000/v1/acquisition-plan \
  -H 'Content-Type: application/json' \
  -d '{
    "atom_groups": [
      {"count": 2, "isotopes": [
        {"shift_mda": 0, "abundance": "3/4"},
        {"shift_mda": 1, "abundance": "1/4"}
      ]}
    ],
    "max_windows": 2,
    "budget": 5
  }'
```

### `POST /v1/detection-probability`

在完整峰型（**不**截取 `top_peaks`）上做**有限离子采样**：`ion_count` 个离子各自独立地
从峰型分布中抽取一个质量偏移；偏移落在窗口 `[l, r]`（闭区间）内即为该窗口命中一次，
**同一离子计入所有包含它的窗口**（窗口可重叠、可相同），落在所有窗口之外的离子
不命中任何窗口（窗外概率保留，不归一化）。某窗口命中数 `≥ min_hits` 即为该窗口**达标**。
返回所有窗口**同时达标**的联合概率，以及按输入顺序排列的各窗口各自达标概率，
全程精确有理数运算，不枚举采样序列、不做随机估计。

请求体（**重叠窗口**示例：偏移 1 同时落入两个窗口）：

```json
{
  "atom_groups": [
    {
      "count": 2,
      "isotopes": [
        {"shift_mda": 0, "abundance": "3/4"},
        {"shift_mda": 1, "abundance": "1/4"}
      ]
    }
  ],
  "ion_count": 2,
  "windows": [
    {"l": 0, "r": 1, "min_hits": 1},
    {"l": 1, "r": 2, "min_hits": 1}
  ]
}
```

字段与约束（任一不满足即返回 **422**）：

| 字段 | 约束 |
| --- | --- |
| `atom_groups` | 与 `/v1/envelope` 完全相同的校验（1～4 组、各组 `count` 总和 ≤ 80、最大可能总偏移 ≤ 2000 等） |
| `ion_count` | 整数，1～80 |
| `windows` | 1～4 个窗口；允许重叠或完全相同 |
| `windows[].l` / `r` | 整数，`0 ≤ l ≤ r ≤ 2000` |
| `windows[].min_hits` | 整数，1～4 |

响应（200）：

```json
{
  "joint_probability": "87/128",
  "marginal_probabilities": ["255/256", "175/256"]
}
```

- `joint_probability`：所有窗口同时达标的概率，已约分的 `"n/d"`（长分数完整输出）。
- `marginal_probabilities`：各窗口各自达标的概率，按**窗口输入顺序**排列，已约分。
- 某窗口不可能达标时（如 `min_hits > ion_count`，或窗口内没有任何峰概率），该窗口的
  达标概率与联合概率精确为 `"0/1"`，其余窗口的边际概率不受影响。

上例峰型为 `0 → 9/16`、`1 → 3/8`、`2 → 1/16`：偏移 1 命中两个窗口，偏移 0 只命中
窗口 0，偏移 2 只命中窗口 1。联合失败只有采样序列 `{0,0}`（81/256）与 `{2,2}`（1/256）
两种，故联合概率为 `1 − 82/256 = 87/128`；窗口 0 仅在 `{2,2}` 时失败（255/256），
窗口 1 仅在 `{0,0}` 时失败（175/256）。

### 调用示例

```bash
curl -s -X POST http://localhost:8000/v1/detection-probability \
  -H 'Content-Type: application/json' \
  -d '{
    "atom_groups": [
      {"count": 2, "isotopes": [
        {"shift_mda": 0, "abundance": "3/4"},
        {"shift_mda": 1, "abundance": "1/4"}
      ]}
    ],
    "ion_count": 2,
    "windows": [
      {"l": 0, "r": 1, "min_hits": 1},
      {"l": 1, "r": 2, "min_hits": 1}
    ]
  }'
```

### `POST /v1/candidate-posterior`

在共同**整数校准偏移**下，对 2～4 个**互斥候选**（各自携带完整峰型与 `n/d` 先验）
做精确贝叶斯推断：先以给定先验抽取一个候选；在该候选下，所有 `ion_count` 个离子
**共用**一个校准偏移 δ（整数，在 `[-calibration_radius, calibration_radius]` 上等先验），
各离子再从该候选的完整峰型中**独立**采样一个原始偏移，其**测得偏移 = 原始偏移 + δ**。
请求给出 1～8 个互异的测得偏移及其离子计数（未列出的偏移全部并入 `other` 箱，
计数为余下离子数）。以观测计数向量的**精确多项似然**对（候选, δ）联合假设做更新，
返回按候选输入顺序、δ 升序排列的联合后验、候选边际后验与偏移边际后验，
全程 `Fraction` 运算，所有分数已约分。

请求体：

```json
{
  "candidates": [
    {
      "atom_groups": [
        {"count": 1, "isotopes": [
          {"shift_mda": 0, "abundance": "1/2"},
          {"shift_mda": 1, "abundance": "1/2"}
        ]}
      ],
      "prior": "1/2"
    },
    {
      "atom_groups": [
        {"count": 1, "isotopes": [
          {"shift_mda": 0, "abundance": "1/4"},
          {"shift_mda": 2, "abundance": "3/4"}
        ]}
      ],
      "prior": "1/2"
    }
  ],
  "ion_count": 2,
  "observations": [
    {"shift_mda": 1, "count": 1},
    {"shift_mda": 2, "count": 1}
  ],
  "calibration_radius": 1
}
```

字段与约束（任一不满足即返回 **422**）：

| 字段 | 约束 |
| --- | --- |
| `candidates` | 2～4 个**互斥**候选（峰型相同也按独立候选处理） |
| `candidates[].atom_groups` | 与 `/v1/envelope` 完全相同的校验（1～4 组、各组 `count` 总和 ≤ 80、最大可能总偏移 ≤ 2000 等） |
| `candidates[].prior` | 与 `abundance` 相同的 `"n/d"` 格式及范围（非负、`≤ 1`）；**各候选先验之和必须精确等于 1** |
| `ion_count` | 整数，1～80 |
| `observations` | 1～8 条，各 `shift_mda` **互异** |
| `observations[].shift_mda` | 整数，0～2000 |
| `observations[].count` | 非负整数，0～80；**各计数之和 ≤ `ion_count`** |
| `calibration_radius` | 整数，0～3；δ 在 `[-r, r]` 上取整数等先验 |

响应（200）。观测在全部假设下概率均为零时，无后验：

```json
{"possible": false}
```

否则：

```json
{
  "possible": true,
  "joint_posteriors": [
    ["0/1", "0/1", "1/1"],
    ["0/1", "0/1", "0/1"]
  ],
  "candidate_posteriors": ["1/1", "0/1"],
  "offset_posteriors": ["0/1", "0/1", "1/1"],
  "offsets": [-1, 0, 1]
}
```

- `offsets`：全部待推断偏移，即 `[-r, r]` 的整数升序列表。
- `joint_posteriors[i][j]`：`P(候选 i, δ=offsets[j] | 观测)`，候选按**输入顺序**、
  组内 δ 按 `offsets` 升序；所有元素之和为 1。
- `candidate_posteriors[i]`：`P(候选 i | 观测)`，候选输入顺序；等于该候选各 δ 列之和。
- `offset_posteriors[j]`：`P(δ=offsets[j] | 观测)`，δ 升序；等于该 δ 各候选之和。
- 多项似然按标准约定 `0^0 = 1`：计数为 0 的箱（包括概率恰为零的箱）不影响似然；
  仅当某个正计数箱（含 `other`）概率为零时该假设似然为零。`possible=false` 表示
  观测在每个（候选, δ）下似然都为零（或全部候选先验为零），此时不返回任何后验。

上例中峰型 A 为 `{0: 1/2, 1: 1/2}`、B 为 `{0: 1/4, 2: 3/4}`，观测到测得偏移
1 与 2 各一个离子：仅 **A + δ=1**（原始偏移 0、1）可能产生，故该联合假设后验为
`1/1`，B 与其余 δ 全部为 `0/1`。

### 调用示例

```bash
curl -s -X POST http://localhost:8000/v1/candidate-posterior \
  -H 'Content-Type: application/json' \
  -d '{
    "candidates": [
      {"atom_groups": [
        {"count": 1, "isotopes": [
          {"shift_mda": 0, "abundance": "1/2"},
          {"shift_mda": 1, "abundance": "1/2"}
        ]}
      ], "prior": "1/2"},
      {"atom_groups": [
        {"count": 1, "isotopes": [
          {"shift_mda": 0, "abundance": "1/4"},
          {"shift_mda": 2, "abundance": "3/4"}
        ]}
      ], "prior": "1/2"}
    ],
    "ion_count": 2,
    "observations": [
      {"shift_mda": 1, "count": 1},
      {"shift_mda": 2, "count": 1}
    ],
    "calibration_radius": 1
  }'
```

### `POST /v1/discrimination-plan`

给定**恰好 2 个等先验候选**（各自由 `atom_groups` 定义完整峰型，**不**截取
`top_peaks`），在 **0～2000 mDa 质量偏移轴**上规划判别窗口：选择至多
`max_windows` 个**互不相交**的整数闭区间 `[l, r]`（**允许空方案**），窗口内每个
格子（含无峰格）都计费，总格数 `total_bins = Σ(r−l+1) ≤ budget`。判别规则：
离子**落窗**时选两候选中**落窗概率较大**者，**未落窗**时选另一方。设两候选的落窗
概率为 `p0`、`p1`（等先验 1/2），判别准确率为 `(1 + |p0 − p1|) / 2`，方案
**精确最大化**该值，全程有理数运算。

并列时依次按以下顺序决胜（与 `/v1/acquisition-plan` 相同）：

1. 准确率相同 → `total_bins` 更少者胜；
2. 仍并列 → 窗口数更少者胜；
3. 仍并列 → 窗口按起点升序、端点展平为 `[l1, r1, l2, r2, …]` 后字典序更小者胜。

请求体：

```json
{
  "candidates": [
    {
      "atom_groups": [
        {"count": 1, "isotopes": [
          {"shift_mda": 0, "abundance": "3/4"},
          {"shift_mda": 1, "abundance": "1/4"}
        ]}
      ]
    },
    {
      "atom_groups": [
        {"count": 1, "isotopes": [
          {"shift_mda": 0, "abundance": "1/4"},
          {"shift_mda": 1, "abundance": "3/4"}
        ]}
      ]
    }
  ],
  "max_windows": 1,
  "budget": 1
}
```

字段与约束（任一不满足即返回 **422**）：

| 字段 | 约束 |
| --- | --- |
| `candidates` | **恰好 2 个**等先验候选（峰型相同也按独立候选处理） |
| `candidates[].atom_groups` | 与 `/v1/envelope` 完全相同的校验（1～4 组、各组 `count` 总和 ≤ 80、最大可能总偏移 ≤ 2000 等） |
| `max_windows` | 整数，1～6；窗口数**至多**为该值 |
| `budget` | 整数，1～200；`total_bins = Σ(r−l+1)` 不得超过该值 |

响应（200）：

```json
{
  "windows": [[0, 0]],
  "total_bins": 1,
  "hit_probabilities": ["3/4", "1/4"],
  "accuracy": "3/4",
  "hit_choice": 0
}
```

- `windows`：`[l, r]` 数组，按起点 `l` 升序；两窗口互不相交；空方案为 `[]`。
- `total_bins`：实际计费格数 `Σ(r−l+1)`，保证 `≤ budget`。
- `hit_probabilities`：两候选各自的落窗概率，按**候选输入顺序**排列，已约分 `"n/d"`。
- `accuracy`：判别准确率 `(1 + |p0 − p1|) / 2`，已约分 `"n/d"`。
- `hit_choice`：落窗时优先的候选索引（`0` 或 `1`）；两候选落窗概率相等时为 `null`。

上例中候选 0 的峰型为 `0 → 3/4, 1 → 1/4`、候选 1 为 `0 → 1/4, 1 → 3/4`：1 格预算下
`[0,0]` 与 `[1,1]` 的 `|p0 − p1|` 都是 1/2（准确率 3/4），字典序决胜选中 `[0,0]`，
落窗时选候选 0。注意目标是**最大化判别准确率**而非覆盖率：上例把 `budget` 放宽到 2
时，`[0,1]` 虽覆盖全部峰（`p0 = p1 = 1`，准确率 1/2），最优方案仍是 `[0,0]`；
两候选峰型完全相同时，任何窗口的 `p0` 都等于 `p1`，空方案（`windows=[]`、
`accuracy="1/2"`、`hit_choice=null`）即为最优。

### 调用示例

```bash
curl -s -X POST http://localhost:8000/v1/discrimination-plan \
  -H 'Content-Type: application/json' \
  -d '{
    "candidates": [
      {"atom_groups": [
        {"count": 1, "isotopes": [
          {"shift_mda": 0, "abundance": "3/4"},
          {"shift_mda": 1, "abundance": "1/4"}
        ]}
      ]},
      {"atom_groups": [
        {"count": 1, "isotopes": [
          {"shift_mda": 0, "abundance": "1/4"},
          {"shift_mda": 1, "abundance": "3/4"}
        ]}
      ]}
    ],
    "max_windows": 1,
    "budget": 1
  }'
```

### `GET /healthz`

健康检查，返回 `{"status": "ok"}`。

## 计算说明

对每组 `count` 个原子，将单原子同位素分布（多项式 `Σ p_i·x^shift_i`）与自身
做 `count` 次精确卷积，得到该组的总偏移分布；再把各组分布依次卷积合并。
所有概率以 `Fraction` 保存，最终约分输出，任意非法输入由 Pydantic 校验并返回 422。

`/v1/mixture-bounds` 以 `w_m = 1 − Σ w_i` 消元后，可行域是权重单纯形与观测区间
截出的有界多面体；通过 `itertools.combinations` 枚举约束子集并做精确高斯消元
得到全部顶点，在顶点上读取各比例的最值与字典序最小可行向量，全程 `Fraction` 运算。

`/v1/peak-attribution` 在逐原子精确卷积之外，为每个部分偏移额外携带各
（组，同位素）的**一阶矩**（概率加权原子数之和）：合并两条路径时，矩按
`m_a·w_b + m_b·w_a` 相加。为避免大分数反复约分，每组丰度分母的最小公倍数
`L_g` 作为公因数提出（该组每条路径带 `L_g**count` 的分母因子），卷积全程在
整数上进行，仅在最后构造一次 `Fraction`：峰概率为 `w / Π L_g**count`，条件期望为
`矩 / w`（公因数约去）；超过目标偏移的部分路径直接舍弃（所有偏移非负）。
目标峰概率为零时返回 `reachable=false`。

`/v1/acquisition-plan` 先复用完整峰型分布（不截取 `top_peaks`），再把所有峰概率
乘以分母的最小公倍数转成**整数**，在峰序列上做动态规划：状态为
（已用窗口数 k、已计费格数 b、最后一个峰是否被开着窗口覆盖），逐峰决定"不覆盖 /
延伸当前窗口（跨过的空格按间隙长度计费）/ 新开一个窗口"。由于把端点滑过无峰格子
只会浪费预算，最优方案的窗口端点必在峰上，因此只需在峰上转移；每个状态内按
（覆盖概率降序、展平端点字典序升序）保留最优前缀——追加公共后缀或同步替换末尾
端点都保持该序，故裁剪是精确的。最终在全部 (k, b) 状态上按
（概率降序、格数升序、窗口数升序、展平端点字典序升序）取最优，覆盖概率以
`Fraction` 一次性约分输出。该带符号权重的窗口 DP 同时被 `/v1/discrimination-plan`
复用（见下）。

`/v1/detection-probability` 先复用完整峰型分布，再把每个偏移映射到其**命中掩码**
（第 j 位置位当且仅当偏移落在窗口 j 内）；沿偏移轴扫描时掩码只在窗口端点处变化，
故单个离子至多有 `2m+1 ≤ 9` 种掩码取值，离子采样即这些掩码上的**多项分布**
（含"不命中任何窗口"的掩码 0，窗外概率因此保留）。随后以各窗口 `min_hits` 为上界
截断的命中计数向量为状态（至多 `5^4 = 625` 个），逐离子做动态规划转移——只追踪
计数向量，**不枚举采样序列**。所有掩码概率乘以其分母的最小公倍数转成整数，DP 全程
在整数上进行（公共分母 `scale**ion_count` 只做符号化携带），联合概率与各窗口边际
概率（从最终状态中按位汇总）各自仅在最后构造一次 `Fraction` 约分输出；不可能达标
的窗口自然得到分子 0，即精确的 `"0/1"`。

`/v1/candidate-posterior` 先为每个候选复用完整峰型分布（不截取 `top_peaks`），
对每个候选与每个 δ ∈ [-r, r]（至多 4 × 7 = 28 个联合假设）把峰型整体平移 δ，
计算观测计数向量的**精确多项概率**：列出的测得偏移各取其平移后峰概率
（缺峰为 0，`0^0 = 1`），未列偏移的总概率 `1 − Σ P(列出偏移)` 作为 `other` 箱
概率、取余下离子数次幂；多项系数对所有假设相同故省略。联合非归一化权重为
`先验_i × (1/(2r+1)) × 多项似然(i, δ)`，直接用 `Fraction` 相乘并各自约分，
求和得到证据 Z；Z = 0 时观测不可能发生，返回 `possible=false`。否则联合后验、
候选边际（按 δ 求和）与偏移边际（按候选求和）均除以同一 Z，各输出一次性约分，
天然归一，不枚举任何离子采样序列。

`/v1/discrimination-plan` 为两个候选各自复用完整峰型分布（不截取 `top_peaks`）。
等先验下最大化准确率 `(1 + |p0 − p1|) / 2` 等价于最大化 `|p0 − p1|`，拆成两个
定向：以 `w0 − w1`（或 `w1 − w0`）为每个偏移格的**带符号权重**，分别套用与
`/v1/acquisition-plan` 相同的窗口 DP——两候选峰概率统一乘以其分母的最小公倍数，
DP 全程在整数上进行；最优窗口端点必落在**权重为正**的格子上（端点滑过非正权重格
不降低捕获权重且严格减少格数），故只需在两峰型支撑集的并集上转移，而"每格保留
（捕获权重降序、展平端点字典序升序）最优前缀"的裁剪对负权重依然精确（追加公共
后缀、同步替换末尾端点都保持该序）。两个定向的最优方案再按
（`|p0 − p1|` 降序、格数升序、窗口数升序、展平端点字典序升序）决胜，定向之间无
偏好；落窗概率、准确率与 `hit_choice` 由最终窗口直接在完整峰型上求和得到，
各分数一次性约分输出。

## 测试

```bash
pip install -r requirements-dev.txt
pytest
```

测试包含与逐路径暴力枚举（`itertools.product`）的对拍、与独立精确求解器
（区间收缩 / 盒子公式 / Fourier-Motzkin 消元）的随机对拍、采集窗口规划与
穷举所有不相交窗口组合的暴力求解器的随机对拍、检出概率与枚举全部离子采样
序列的暴力求解器的随机对拍、候选后验与枚举每个（候选, δ）下全部离子序列的
贝叶斯暴力求解器的随机对拍、判别窗口规划与穷举所有不相交窗口组合并按
（|p0−p1|、格数、窗口数、展平端点）决胜的暴力求解器的随机对拍，
以及全部非法输入的 422 用例。

## 项目结构

```
app/main.py          # FastAPI 应用：校验模型 + 精确卷积 + 精确有理数 LP + 带符号权重窗口 DP（采集/判别共用）+ 检出概率 DP + 候选/偏移后验
tests/test_envelope.py
tests/test_peak_attribution.py
tests/test_mixture_bounds.py
tests/test_acquisition_plan.py
tests/test_detection_probability.py
tests/test_candidate_posterior.py
tests/test_discrimination_plan.py
requirements.txt     # 运行依赖
requirements-dev.txt # 测试依赖
Dockerfile
docker-compose.yml
```
