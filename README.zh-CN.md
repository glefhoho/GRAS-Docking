# GRAS-Docking（中文说明）

把食品添加剂（GRAS 名单化学品）批量对接到一组蛋白结构上，并按**靶点特异性**
而不是原始对接分数来排序的四步流程。完整说明见英文
[README.md](README.md)。

> **重要前提**：这个流程产出的是"值得做实验验证的假设"，不是结合证据。
> 在解读任何结果之前，请先读 [`docs/CAVEATS.md`](docs/CAVEATS.md)。

## 四步流程

| 步骤 | 脚本 | 输入 → 输出 |
|---|---|---|
| 1 | `scripts/step1_fetch_smiles.py` | `Name,CAS` 的 CSV → `chemicals_prepared.csv`（PubChem 查询、脱盐、校验） |
| 2 | `scripts/step2_build_ligands.py` | SMILES → `ligands/*.pdbqt`（Open Babel 生成 3D 构象并按 pH 加质子） |
| 3 | `scripts/step3_dock.py` | 蛋白 + 配体 → `docking_results.csv`（smina 对接，口袋来自 fpocket） |
| 4 | `scripts/step4_analyze.py` | 打分矩阵 → `druggable.csv`（双向 z 分数 + 配体效率） |

**排序思路**：原始对接分数主要由分子大小决定，直接横向比较等于按分子量排序。
所以每个分数做两次标准化——对该化合物自身在所有蛋白上的均值（`Z_Chem`：
它是不是特别"偏爱"这个靶点），以及对该蛋白在所有化合物上的均值（`Z_Prot`：
这个蛋白是不是对谁都亲和）。真正有意义的特异性结合应该在两个方向上都突出。

## 安装

```bash
conda env create -f environment.yml
conda activate gras-docking
```

三个核心程序（`obabel` / `smina` / `fpocket`）不能用 pip 装，必须走 conda 或
自行编译。检查：

```bash
obabel -V && smina --version && fpocket -h | head -1
```

RDKit 不是必需，但强烈建议装——没有它，脱盐会退化成较弱的字符串启发式，
而且算不了配体效率。

## 快速开始

```bash
cp data/example_chemicals.csv chemicals_list.csv

python scripts/step1_fetch_smiles.py -i chemicals_list.csv
python scripts/step2_build_ligands.py --ph 6.5

mkdir -p proteins                          # 把 .pdb 放进来
for f in proteins/*.pdb; do fpocket -f "$f"; done

python scripts/step3_dock.py --exhaustiveness 16 --cpu 8
python scripts/step4_analyze.py
```

前三步都支持断点续跑：中断后重新执行同一条命令，已完成的部分会跳过；
加 `--overwrite` 才会重头再来。

## 几个必须自己确认的参数

* `--ph`（步骤 2）：按你研究体系的实际 pH 设置，默认 6.5。苯甲酸这类
  羧酸在不同 pH 下的解离状态完全不同，对接结果也不同。
* `--pocket-rank`（步骤 3）：fpocket 排第一的口袋**不一定**是活性位点。
  每个蛋白都应该对照已知催化残基或共晶配体确认一次。
* `--exhaustiveness`（步骤 3）：默认 8 只适合试跑，要发表请用 16–32。
* `--score-lead-min`（步骤 4）：默认 −8 kcal/mol 这个绝对阈值**带有明显的
  分子大小偏倚**，小分子（如苯甲酸）几乎不可能达到。可以用
  `--no-lead-filter` 关掉，只按特异性排序。

## 目录结构

```
GRAS-Docking/
├── scripts/          维护中的正式脚本
├── legacy/           原始四个脚本，一字未改，用于溯源
├── data/             示例输入
├── docs/
│   ├── CHANGES.md    相对原脚本改了什么、哪些改动会影响结果
│   └── CAVEATS.md    方法学局限，发表前必读
├── tests/            pytest 测试（不需要联网和外部程序）
├── environment.yml
└── requirements.txt
```

`pytest -q` 可以直接跑测试。

## 许可

MIT，见 [LICENSE](LICENSE)。
