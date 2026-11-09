<p align="center"><img src="docs/logo.png" alt="AnnoLens" width="160"></p>

<h1 align="center">AnnoLens</h1>

<p align="center">像人工审阅员一样看证据、修基因模型——智能体 + 技能，以 MCP 服务提供，适用于任何物种</p>

<p align="center"><b>中文</b> | <a href="README.en.md">English</a></p>

---

## 简介

自动注释流程（EVM、BRAKER、MAKER 等）给出的基因模型里，总有一部分有问题：截短、被拆成两段、两个基因被并成一个、缺起始或终止密码子、没有 UTR；还有一些基因根本没被注释出来。人工审阅员会在基因组浏览器里对照转录本、同源蛋白和剪接读段，逐个改正。AnnoLens 把这个过程交给能看图的大模型来做，同时保证它不能乱改：

- **智能体看证据图做判断**：为每个基因或位点画一张证据图，内容包括当前模型、各样本的观测转录本、同源蛋白在基因组上的投影和剪接读段，并附一份文字证据表。
- **技能规定怎么审**：先判断基因的范围，再确定结构；哪些证据可信，什么时候保留、修改、新增或判为无法判断，都写在版本化的技能文本里。
- **引擎守住硬规则**：候选模型只能取自观测转录本或完整的同源蛋白投影，模型不能自己写坐标。完整 ORF、经典剪接位点、不与相邻基因冲突、内含子有读段或投影支持、长度与同源蛋白相符、不能比多数样本支持的 ORF 更短，这些都由引擎检查，模型说了不算。
- **只提建议，不改原件**：输出是审阅记录和候选 GFF3，输入的注释和基因组都不会被修改。

## 功能

| 功能 | 说明 |
|---|---|
| 全基因组规则初查 | `annolens scan` 逐个检查全部基因：编码完整性、UTR、基因间重叠、证据支持、与独立注释比长短、拆分与融合、读段证实的移码、重复序列重叠，并找出没有模型的位点。有问题的基因按优先级排进队列 |
| 修复已有基因 | 截短的延长；被拆开的合并；缺起始或终止的补全；UTR 取自观测转录本 |
| 新增缺失基因 | 位点有转录本或完整同源投影、却没有模型时，从证据构建新基因 |
| 读段核查 | 在任意位置统计基因组长读长（如 HiFi）和 RNA 读段是否支持组装序列，以区分组装错误和真实存在的移码或终止 |
| 规则补 UTR | 用内含子链一致的观测转录本批量补 UTR，并检查是否与相邻基因冲突 |
| 批量审查 | 并行审几千到几万个基因；可接任何兼容 OpenAI 或 Anthropic 接口、能看图的模型 |
| MCP 服务 | 供 Codex、Claude Code、WorkBuddy、Cursor 等支持 MCP 的宿主调用 |

## 安装

服务器端（存放数据的机器）：

```bash
pip install git+https://github.com/forageseed/AnnoLens.git
# 外部程序：miniprot（同源投影）；只有 BAM、没有 GTF 时还需要 StringTie
```

## 快速开始

1. 复制 [`examples/annolens.toml`](examples/annolens.toml)，填入你的基因组、注释、RNA 比对、同源蛋白（可选：长读长转录本、HiFi、重复序列注释、模型接口）。
2. 生成工作目录，并做全基因组初查：

```bash
annolens prepare --manifest annolens.toml --threads 16
annolens scan    --manifest annolens.toml
```

3. 批量审查（模型密钥只从环境变量读取）：

```bash
export ANNOLENS_API_KEY=...        # 变量名由清单里的 api_key_env 决定
annolens review  --manifest annolens.toml --queue P1 --workers 4
annolens addgene --manifest annolens.toml --queue missing_high --workers 4
```

结果写在 `<work_dir>/reviews/`：每个基因一条审阅记录（JSON），通过验证的修改和新增基因同时输出为 GFF3。

### 同源蛋白的三种角色

| role | 含义 | 用法 |
|---|---|---|
| `independent` | 同一物种另一套独立产生的注释 | 锚定起点、终点和长度；判断模型是否截短的主要参照 |
| `related` | 近缘种的蛋白组 | 锚定起点、终点和长度 |
| `paralog` | 与被审注释出自同一流程的蛋白，例如同源拷贝 | 只用来支持内含子 |

不要把“当初用作注释证据的蛋白组”设为 `independent`：它会重复注释里同样的截短。

## 在 Codex、Claude Code、WorkBuddy 等宿主中使用（MCP）

一条命令就能同时**登记 MCP 服务并安装技能**，两者来自同一份技能文本，不会对不上：

```bash
annolens install --host codex  --ssh <服务器别名> --manifest /path/on/server/annolens.toml
annolens install --host claude --ssh <服务器别名> --manifest /path/on/server/annolens.toml
annolens install --host generic --ssh <服务器别名> --manifest /path/on/server/annolens.toml   # 打印 mcpServers JSON
```

数据在本机时可以省略 `--ssh`。改动宿主配置之前，安装器会先备份为 `*.annolens.bak`。Windows 用户请在 PowerShell 中运行。

装好后直接用自然语言提问即可：

- “帮我看看这个基因模型对不对：<gene_id>”
- “这个位置是不是漏注释了一个基因？<chrom>:<start>-<end>”
- “把 P1 队列的前 20 个基因审一遍”

| 步骤 | 已有基因 | 无模型位点 |
|---|---|---|
| 列队列 | `annolens_list_queue` | `annolens_list_queue` |
| 打开 | `annolens_open_gene` | `annolens_open_locus` |
| 换窗口 | `annolens_gene_view` | `annolens_locus_view` |
| 列 ORF | `annolens_gene_orfs` | `annolens_locus_orfs` |
| 构建候选 | `annolens_gene_build` | `annolens_locus_build`（用转录本或同源投影） |
| 读段核查 | `annolens_readcheck` | `annolens_readcheck` |
| 下结论 | `annolens_gene_decide` | `annolens_locus_decide` |

另外提供两个提示词 `annolens_review_gene`、`annolens_add_gene`，以及资源 `annolens://skill/gene`、`annolens://skill/locus`。每个会话能查看、构建、核查读段几次，由服务端限定；判为修改或新增之前，服务端会再验证一遍。

## 结论类型

| 已有基因 | 含义 |
|---|---|
| `retain` | 当前模型完整，且有证据支持 |
| `recommend_candidate` | 采用经引擎验证的新结构 |
| `genome_frameshift` | 独立蛋白需要的移码或终止，被读段证实确实存在于本基因组；保留原模型并做标记 |
| `assembly_conflict` | 组装层面的问题（缺口、错接，或读段显示组装有误） |
| `unresolved` | 证据冲突，无法判断 |

| 无模型位点 | 含义 |
|---|---|
| `add_gene` | 新增经引擎验证的基因 |
| `no_gene` | 不是蛋白编码基因（转座子、非编码转录本或证据不足） |
| `genome_frameshift` | 该基因在本基因组中已经失活 |
| `unresolved` | 无法判断 |

## 设计要点

- **结构只能取自观测证据**：模型只能选“用哪条转录本”“哪个 ATG”“用哪条同源投影”，不能直接写坐标。
- **多证据门控**：有长读长全链支持；或者内含子有读段支持、起止和长度有同源投影锚定；或者是多数样本一致的主导 ORF。至少满足其中一条，修改才能被采纳。
- **只修确实有问题的模型**：完整且已有证据支持的模型不会被替换。唯一的例外是多样本主导 ORF 和独立注释的完整蛋白同时表明它被截短了。
- **防截短**：长度比主导 ORF 短、支持样本也更少的候选，一律不能采纳。
- **可追溯**：每个候选都带内容哈希、引擎哈希和输入指纹；证据或代码一变，旧候选自动作废。

## 已知限制

- 需要宿主模型能看图；只有文字证据表时，判断能力会下降。
- 目前还不能把一个基因拆成两个。
- 规则补 UTR 依赖转录本内含子链与模型完全一致；拿不准的情况会交给智能体处理。

## 目录

```
annolens/            Python 包：config、engine、tools、addgene、readcheck、scan、prepare、session、agent、mcp_server、install、render
annolens/skills/     版本化的技能文本（已有基因 gene_v3.5、新增基因 locus_v3.6）
examples/            项目清单示例
scripts/             发布前脱敏检查 check_sensitive.py
```

## 引用与许可

论文准备中。许可证：[Apache-2.0](LICENSE)。
