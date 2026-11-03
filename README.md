# CBF 提取（Native Space）

在**个体原生空间（native space）**批量提取 ASL-CBF 的脑区平均值，同时使用 **AAL3** 和 **LPBA40** 两套模板。

与常规做法（先把 CBF 标准化到 MNI 再提 ROI）不同，本脚本**反过来把图谱变形到每个被试自己的 T1 空间**，因此不需要对 CBF 图做重采样平滑，避免了标准化插值带来的部分容积效应。

## 处理流程

对每一次扫描：

1. **配准** — 用 `spm_coreg`（归一化互信息）把 CBF 图配准到同一被试的 T1
2. **逆变形** — 用分割得到的形变场 `y_T1.nii`，把 AAL3 和 LPBA40 两张图谱从 MNI 拉回该被试的 T1 空间
3. **重采样** — 把图谱和灰质概率图重采样到 CBF 的体素网格（图谱用最近邻保持标签整数，灰质用线性插值）
4. **提取** — 在「属于该 ROI **且** 灰质概率 ≥ 0.30」的体素内取 CBF 平均值
5. **存盘** — 每处理完一次扫描就写一次 checkpoint

## 环境要求

- **MATLAB**（较新版本，脚本用了 `string`、`erase` 等）
- **SPM**（脚本头部标注在 SPM25 下运行；SPM12 的 `spm.util.defs` 批处理接口同样适用）

模板文件本身见 [brainatlases-](https://github.com/weishenshuying/brainatlases-) 仓库。

## 使用前必须修改的地方

脚本里的路径是硬编码的，**换机器必须改**。集中在开头：

| 行 | 变量 | 改成什么 |
|---|---|---|
| 20 | `addpath('/Users/mac/spm')` | 你的 SPM 安装路径 |
| 27 | `aslDir = '/Users/mac/EXP/eTRE/ASL'` | ASL 数据的根目录 |
| 30–31 | `lpba` | LPBA40 标签图路径（用的是 `maxprob/` 里的灰质标注图） |
| 33 | `aal` | AAL3 模板路径（1 mm 版本 `AAL3v1_1mm.nii`） |

> `runDir` 是 `aslDir` 下的 `asl_run`，扫描目录都放在这一层。

如果你的数据集不是这批，**还要改样本相关的参数**（脚本里有断言会拦住你）：

| 位置 | 含义 | 当前值 |
|---|---|---|
| 69 | `excludePatients` | 排除 `002, 010, 053, 060` |
| 71–72 | `excludeH` | 排除 `H014, H018, H021, H032, H052, H072` |
| 64 | M0 / M1 目录数 | 各 44 个 |
| 79 / 82 | 入组患者 / 对照数 | 40 / 42 |

## 输入目录结构

```
<aslDir>/asl_run/
├── 010_M0/
│   ├── STRUC/
│   │   ├── T1.nii
│   │   └── segmentation_work/
│   │       ├── c1T1.nii          ← 灰质概率图（SPM 分割产物）
│   │       └── y_T1.nii          ← 到 MNI 的形变场
│   └── FUNC/
│       └── meanCBF_0_sASLflt_rPASL.nii    ← CBF 图
├── 010_M1/          （同结构）
├── 011_M0/ …
└── H001/            （健康对照，同结构）
```

命名规则：患者目录是 `3位编号_M0` / `3位编号_M1`，对照目录是 `H` + 3 位数字。

运行前脚本会**先检查全部 122×4 个文件是否齐全**，缺任何一个都会列出清单并中止，不会跑一半才报错。

## 输出

在 `<aslDir>/NativeCBF_122Scans_<时间戳>/` 下：

| 文件 | 内容 |
|---|---|
| `LPBA40_Native_CBF_40Pairs_42H.xlsx` | 122 行 × 55 列，LPBA40 各 ROI 的平均 CBF |
| `AAL3_Native_CBF_40Pairs_42H.xlsx` | 122 行 × 150 列，AAL3 各 ROI 的平均 CBF |
| `Processing_Log.csv` | 每次扫描的成功与否、有效 ROI 数、配准平移旋转量、错误信息 |
| `Native_CBF_Results.mat` | 全部变量的 MAT 存档 |
| `Progress_Checkpoint.mat` | 每扫一次刷新一次，中断后可续跑 |
| `<扫描ID>/` | 该扫描的中间产物：`Native_Coreg_CBF.nii`、逆变形后的图谱、`*_on_CBF.nii`（图谱与灰质在 CBF 网格上的版本，可直接查看核对） |

两张 Excel 的前三列是 `Group`（M0/M1/H）、`Patient_ID`、`Scan_ID`，后面是 `ROI_<编号>`，编号即图谱原始的 ROI 整数值（不是连续重排的 1–N），方便直接对照字典。

## 关键参数

| 参数 | 值 | 位置 |
|---|---|---|
| 灰质阈值 | `G >= 0.30` | 390 / 412 行 |
| LPBA40 ROI 数 | 55 个（局部 ID 列表见 154–155 行） | 154 |
| AAL3 ROI 数 | 150 个（局部 ID 列表见 157–158 行） | 157 |
| 配准代价函数 | 归一化互信息，`sep=[4 2]`，`fwhm=[7 7]` | 194–198 |
| 图谱重采样 | 最近邻（保持标签为整数） | 318–321 |
| 灰质重采样 | 线性插值，后截断到 [0,1] | 同上 |

## 注意事项

- **复核标志**：平移超过 30 mm 或旋转超过 20° 的扫描会在 `Processing_Log.csv` 里标 `NeedsReview=true`，这些多半是配准失败，建议人工看一眼再决定是否纳入统计。
- **原始文件不会被修改**：CBF 图是先 `copyfile` 到输出目录再改写空间信息的。
- **脚本开头有 `clear; clc;`**，会清空工作区，别在保留变量的情况下直接跑。
- **ROI 编号是图谱原值**：AAL3 的 150 个是从 166 个区域里筛出来的，不是 1–150 连续编号，做统计时别按列序号去查字典。

## 引用

使用本脚本产出的结果，请引用原始图谱文献：

- **AAL3** — Rolls, E. T., Huang, C. C., Lin, C. P., Feng, J., & Joliot, M. (2020). Automated anatomical labelling atlas 3. *NeuroImage*, 206, 116189.
- **LPBA40** — Shattuck, D. W., Mirza, M., Adisetiyo, V., et al. (2008). Construction of a 3D probabilistic atlas of human cortical structures. *NeuroImage*, 39(3), 1064-1080.
