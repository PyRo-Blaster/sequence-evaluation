# Sequence Evaluation

`sequence-evaluation` 是一个面向蛋白序列分析的 agent skill，重点覆盖抗体和治疗性蛋白的常见评估任务。

它适合在用户提供 FASTA 文件或氨基酸序列后，快速完成以下工作：

- 计算分子量（MW）、等电点（pI）、消光系数 / A280 以及可开发性指标
- 识别抗体可变区 CDR（优先使用 abnumber/ANARCI，回退到锚点正则）
- 映射 IgG Fc 突变并按 EU 编号解释，自动标注常见工程突变
- 扫描可开发性序列风险（糖基化位点、脱酰胺、异构化、氧化等）
- 在已有预测结构上标记 charged patch 与疏水表面 patch（不是 TAP / SAP）
- 对变量区或一般蛋白做同源序列搜索
- 组织抗体-抗原复合物结构预测
- 分析 paratope / epitope 接触残基，并映射到 CDR
- 由各阶段 JSON 输出确定性地汇总为结构化 Markdown 报告

## 适用场景

这个 skill 尤其适合以下分子类型：

- 单抗、双抗、三抗
- VHH / nanobody
- 含 Fc 工程改造的治疗性抗体
- 需要基础理化分析和同源搜索的一般蛋白

典型触发问题包括：

- `analyze this FASTA`
- `evaluate my antibody sequences`
- `find the CDRs`
- `what Fc mutations are present`
- `run a homology search on my protein`

## 仓库内容

```text
.
├── SKILL.md
├── README.md
├── LICENSE
├── examples/
│   ├── trastuzumab.fasta        # 测试用曲妥珠单抗序列
│   ├── mini_complex.cif         # 接触分析用小型复合物
│   └── mini_cdrs.json
├── references/
│   ├── antibody-numbering.md
│   ├── developability.md
│   └── structural-prediction.md
├── tests/
│   └── smoke_test.sh            # 端到端冒烟测试
└── scripts/
    ├── _common.py               # 共用工具（序列校验 / JSON）
    ├── analyze_interfaces.py
    ├── analyze_properties.py
    ├── build_cdr_map.py         # 由 find_cdrs 输出拼装 CDR map（三处共用）
    ├── select_structure_model.py # 按模态选结构预测模型（VHH→Boltz-1 / Fab→Chai-1）
    ├── compile_report.py        # 由 JSON 汇总分阶段 Markdown 报告
    ├── find_cdrs.py
    ├── find_mutations.py
    ├── scan_liabilities.py      # 可开发性序列风险扫描（支持 CDR 升级）
    └── scan_surface_patches.py  # 结构表面 charged / hydrophobic patch
```

所有脚本均支持 `--json` 输出，便于用 `compile_report.py` 自动汇总报告。
运行 `bash tests/smoke_test.sh` 可在新环境中验证脚本和依赖是否可用。

## Skill 能力概览

### 1. 理化性质分析

使用 `scripts/analyze_properties.py` 计算：

- 单链 MW / pI
- 消光系数与 A280（1 g/L，用于浓度测定）
- 不稳定指数、脂肪族指数、GRAVY 等可开发性指标
- 去除信号肽后的成熟链性质
- 多链复合体总 MW / pI
- 二硫键对复合体分子量的修正
- 对非标准氨基酸（X/B/Z/U…）做输入校验，避免计算中途崩溃

示例：

```bash
uv run scripts/analyze_properties.py input.fasta --complex --disulfide-bonds 18
```

### 2. 抗体区段注释

使用 `scripts/find_cdrs.py` 识别：

- VH CDR1 / CDR2 / CDR3
- VL CDR1 / CDR2 / CDR3
- VHH 的 CDR 区段

示例：

```bash
uv run scripts/find_cdrs.py --vh EVQLVES... --name VH_A
uv run scripts/find_cdrs.py --vl DIQMTQ... --name VL_A
```

### 3. Fc 突变映射

使用 `scripts/find_mutations.py`：

- 自动检测同种型（IgG1 / IgG2 / IgG4），与对应野生型比对
- 使用 EU 编号输出突变位点（基于全局比对，可正确处理插入/缺失）
- 自动标注常见工程位点，如 `LALA-PG`、`LS`、`YTE`、`KiH`、IgG4 `S228P`
- 可用 `--isotype` 指定同种型，或用 `--reference-fasta` 提供权威参考序列

示例：

```bash
uv run scripts/find_mutations.py "MGWSCIILFLV...ASTKGPSVF..." --label HC1
uv run scripts/find_mutations.py "<seq>" --isotype IgG4
```

> 说明：内置的 IgG2/IgG4 参考序列依据同种型典型差异重建；仅 IgG1 经审批药物验证。
> 用于正式/法规场景时建议通过 `--reference-fasta` 提供 IMGT/UniProt 权威序列。

### 4. 可开发性风险扫描

使用 `scripts/scan_liabilities.py` 标注序列层面的化学/翻译后修饰风险：

- N-糖基化位点（N-X-S/T）
- 脱酰胺（NG/NS…）、异构化（DG/DS…）
- Asp-Pro 断裂、Met/Trp 氧化热点
- 未配对半胱氨酸、N 端焦谷氨酸

示例：

```bash
uv run scripts/build_cdr_map.py HC=vh.json LC=vl.json > cdrmap.json
uv run scripts/scan_liabilities.py input.fasta --cdr-json cdrmap.json --json
```

传 `--cdr-json` 后，落在 CDR 内的风险自动升一级（氧化位点若还落在有接触的 CDR
上再升一级）；不传则全部是默认等级，并会在输出里说明。

有预测结构之后，再用 `scripts/scan_surface_patches.py` 扫抗体链（不要包含抗原）上的正/负电荷 patch 和疏水表面 patch。这是 `surface-patch-v1`，不是 TAP，也不是 SAP。没有结构时不要从序列推断 patch。

```bash
uv run scripts/scan_surface_patches.py pred.model_idx_0.cif \
  --chains A B --cdr-json cdrs.json --json
```

### 5. 同源搜索

结合 `protein-sequence-similarity-search` skill：

- 对 VH / VL / VHH 单独做搜索
- 汇总 top hits、覆盖度、E-value 和 identity
- 辅助推断来源、家族和相似已知分子

### 6. 结构预测与界面分析

结构预测模型按模态固定选择：**VHH / nanobody → Boltz-1**（`boltz` skill），
**Fab / Fv（VH+VL）→ Chai-1**（`chai` skill）。用脚本解析并准备 FASTA，
同时固定链顺序（A = 第一条记录）：

```bash
uv run scripts/select_structure_model.py --cdrs vh.json vl.json \
  --seqs VH=vh.seq VL=vl.seq --out-fasta fv.fasta
```

默认预测 **apo** Fv/VHH：表面 patch 关心的是分子在制剂里暴露的表面，
不是被抗原埋住的界面。界面分析（paratope/epitope）不再是报告阶段，
需要抗原复合物时单独运行：

结合 `chai` 和 `modal` skill：

- 为抗体-抗原复合物准备多链 FASTA
- 运行 Chai-1 预测复合物结构
- 用 `scripts/analyze_interfaces.py` 解析 paratope / epitope 接触

示例：

```bash
uv run scripts/analyze_interfaces.py pred.model_idx_0.cif \
  --antibody A B \
  --antigen C \
  --cutoff 5.0
```

## 环境预检与后端选择

每个会话开始时先运行一次环境预检，它会一次性检查**工具链、Python 包和 GPU 算力**：

```bash
uv run scripts/preflight.py                  # 快速检查
uv run scripts/preflight.py --provision      # 另外验证 uv 能拉取核心依赖
uv run scripts/preflight.py --check-network  # 另外探测外网连通性
```

检查内容：`uv` 是否可用；Modal CLI 是否安装**并已认证**；`biopython`/`numpy`（由 uv
按脚本内联依赖自动安装）、可选的 `abnumber`、以及 `modal`/`torch`/`chai_lab`/`boltz`；
GPU 型号与显存、CPU 核数、内存。缺什么会直接给出安装命令；所有必需项满足时退出码为 0。

它会报告 GPU（NVIDIA / Apple MPS）、CPU、内存、Modal CLI 与网络状况，并按任务类型给出建议：

- 轻量 CPU 任务（理化性质、CDR、突变、可开发性、界面几何）→ **始终本地**
- 抗体结构建模（IgFold/ImmuneBuilder 等）→ 有 GPU 用本地；否则 CPU（较慢）或远程
- 复合物 cofolding（Chai-1/Boltz/AF3，**需要 GPU**）→ 本地显存足够则本地；否则 **Modal**；再否则 **Web API/服务器**

### Modal 云端 GPU（Chai-1 / Boltz）

```bash
pip install modal && modal token new      # 一次性设置
modal run modal/modal_chai1.py --input-faa in.fasta --out-dir predictions/x
modal run modal/modal_boltz.py --input-path in.fasta --out-dir predictions/y
```

用 `MODAL_GPU` 选择 GPU（默认 `A100-40GB`）。权重缓存在持久化 Modal Volume 中。
Chai-1 不使用 MSA 服务器；Boltz 的 `--use-msa-server` 默认关闭（会将序列发送到公共
ColabFold 服务器）。详见 [modal-setup.md](./references/modal-setup.md)。

> 原则：GPU 密集任务在本机不具备条件时，使用 Modal 或在线服务（如 AlphaFold3 server）；
> 对敏感序列优先本地或 Modal，避免发送到第三方。

## 推荐工作流

对于抗体序列，通常按这个顺序使用：

1. 运行理化性质分析（多聚体务必传 `--stoichiometry`）
2. 识别 VH / VL / VHH 的 CDR
3. 用 `build_cdr_map.py` 拼装 CDR map（后面三个脚本共用）
4. 检查 Fc 突变及功能含义
5. 扫描可开发性序列风险，并传 `--cdr-json` 让 CDR 内的风险自动升级
6. 对变量区做同源搜索
7. 如需结构层面解释，再做复合物预测、界面分析和表面 patch 扫描
8. 用 `compile_report.py` 汇总为分阶段报告

## 分阶段报告

`compile_report.py` 输出 9 个阶段，**每个阶段都给出数据（Findings）和结论
（Assessment）**，没有输入的阶段会显式写明 `Not assessed`，而不是悄悄省略：

| 阶段 | 结论回答什么 |
|---|---|
| 1 Molecule Triage | 这是什么分子，链组成是否自洽 |
| 2 Chain Properties | 有没有超出配方/表达的常规范围 |
| 3 Domain Annotation | CDR 边界可不可信，用的哪个引擎 |
| 4 Engineering Audit | 突变集做什么用；KiH 是否互补 |
| 5 Sequence Liabilities | 哪些风险落在 CDR 内，先改哪个 |
| 6 Structural Developability | 哪些 patch 碰到 CDR，对模型的依赖有多大 |
| 7 Lineage | 来源与最接近的已知分子 |
| 8 Verdict | 风险清单排序 + STOP / ENGINEER / PROCEED |

第 8 阶段汇总前 7 阶段提出的所有风险，并给出明确结论；发现构建体错误
（例如 knob 没有对应的 hole）时直接给 **STOP**。报告末尾的 Provenance 记录
时间戳、各输入文件、复合体化学计量、patch 参数、CDR 引擎与 Fc 参考同种型。

> Assessment 层是由 `references/` 里的阈值规则推导出来的，不要改成自由发挥的
> 叙述，也不要加 JSON 支持不了的结论。

对于一般蛋白，通常只需要：

1. 理化性质分析
2. 同源搜索

## 参考文档

- [SKILL.md](./SKILL.md): skill 主说明
- [antibody-numbering.md](./references/antibody-numbering.md): Kabat / IMGT / EU 编号说明
- [developability.md](./references/developability.md): 可开发性风险说明
- [modal-setup.md](./references/modal-setup.md): Modal 云端 GPU 配置（Chai-1 / Boltz）
- [structural-prediction.md](./references/structural-prediction.md): Chai-1 结构预测与界面分析说明

## 使用建议

- 抗体相关输入尽量区分 `VH`、`VL`、`VHH` 和完整 heavy chain
- 做 CDR 分析时优先传入可变区，不要混入常数区
- 做 Fc 突变分析时应传完整 heavy chain，便于自动定位常数区锚点
- 做复合体预测时，链顺序要和后续界面分析参数保持一致

## License

本项目以 [MIT License](./LICENSE) 发布。
