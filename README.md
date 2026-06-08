# Sequence Evaluation

`sequence-evaluation` 是一个面向蛋白序列分析的 agent skill，重点覆盖抗体和治疗性蛋白的常见评估任务。

它适合在用户提供 FASTA 文件或氨基酸序列后，快速完成以下工作：

- 计算分子量（MW）和等电点（pI）
- 识别抗体可变区 CDR
- 映射 IgG Fc 突变并按 EU 编号解释
- 对变量区或一般蛋白做同源序列搜索
- 组织抗体-抗原复合物结构预测
- 分析 paratope / epitope 接触残基
- 汇总为结构化 Markdown 报告

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
├── references/
│   ├── antibody-numbering.md
│   └── structural-prediction.md
└── scripts/
    ├── analyze_interfaces.py
    ├── analyze_properties.py
    ├── find_cdrs.py
    └── find_mutations.py
```

## Skill 能力概览

### 1. 理化性质分析

使用 `scripts/analyze_properties.py` 计算：

- 单链 MW / pI
- 去除信号肽后的成熟链性质
- 多链复合体总 MW / pI
- 二硫键对复合体分子量的修正

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

- 对照 WT `IGHG1*01`
- 使用 EU 编号输出突变位点
- 方便标注常见工程位点，如 `LALA-PG`、`LS`、`YTE`、`KiH`

示例：

```bash
uv run scripts/find_mutations.py "MGWSCIILFLV...ASTKGPSVF..." --label HC1
```

### 4. 同源搜索

结合 `protein-sequence-similarity-search` skill：

- 对 VH / VL / VHH 单独做搜索
- 汇总 top hits、覆盖度、E-value 和 identity
- 辅助推断来源、家族和相似已知分子

### 5. 结构预测与界面分析

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

## 推荐工作流

对于抗体序列，通常按这个顺序使用：

1. 运行理化性质分析
2. 识别 VH / VL / VHH 的 CDR
3. 检查 Fc 突变及功能含义
4. 对变量区做同源搜索
5. 如需结构层面解释，再做复合物预测和界面分析
6. 输出完整 Markdown 报告

对于一般蛋白，通常只需要：

1. 理化性质分析
2. 同源搜索

## 参考文档

- [SKILL.md](./SKILL.md): skill 主说明
- [antibody-numbering.md](./references/antibody-numbering.md): Kabat / IMGT / EU 编号说明
- [structural-prediction.md](./references/structural-prediction.md): Chai-1 结构预测与界面分析说明

## 使用建议

- 抗体相关输入尽量区分 `VH`、`VL`、`VHH` 和完整 heavy chain
- 做 CDR 分析时优先传入可变区，不要混入常数区
- 做 Fc 突变分析时应传完整 heavy chain，便于自动定位常数区锚点
- 做复合体预测时，链顺序要和后续界面分析参数保持一致

## License

如果你准备将这个 skill 进一步公开分发，建议补充一个明确的开源许可证，例如 `MIT`。
