# fusion-metrics（中文）

[English](README.md) · [Tiếng Việt](README.vi.md)

**图像融合客观评价指标的 Python 实现：从 MATLAB 逐行移植，并与原版逐值核对到舍入误差级别。**

作者：**[Vu Viet Anh（武越英）](https://vuvietanh2309.github.io/)**，AI 工程师，越南河内百科大学（HUST）计算机科学硕士在读

## 为什么需要这个库

图像融合论文通常报告 EN、SD、MI、VIF、Q<sup>AB/F</sup>、SSIM 等指标，但每个团队使用的 MATLAB 脚本或 Python 改写版都不一样，
细节上的差异往往没人注意：灰度转换系数不同、`conv2(...,'same')` 错位一个像素、`uint8` 饱和截断、`round` 规则、`edge()`
细化方式……结果是同一个方法在不同论文里得到不同的数值，表格之间无法直接比较。

`fusion-metrics` 提供**一个统一、可 `pip install` 的实现**：

- **忠实于 MATLAB。** 每个函数都从一个参考 `.m` 文件逐行移植，保留函数名、参数顺序和数值细节。
- **经过验证。** 自检脚本把每个值与原版 `.m` 文件（在 GNU Octave 中运行）的结果逐一比较，相对误差 ≤ 1e-12；
  CI 在 Python 3.10–3.13 上自动运行。
- **可直接用于基准测试。** 既可以在 Python 中计算单张图像，也可以用命令行批量评估一个方法的全部结果并导出 CSV。

## 安装与快速上手

```bash
pip install fusion-metrics
```

命令行（按文件名配对 `ir/`、`vi/` 与 `fused/` 中的图像）：

```bash
fusion-metrics --dirs ir/ vi/ fused/ --preset ivif --csv results.csv
fusion-metrics --list          # 查看全部指标和预设
```

Python：

```python
import fusion_metrics as fm

A = fm.read_gray("ir.png")      # 与 MATLAB imread + rgb2gray 完全一致
B = fm.read_gray("vi.png")
F = fm.read_gray("fused.png")

print(fm.evaluate(A, B, F, "ivif"))
```

## 指标与预设

| 预设 | 指标 |
|---|---|
| `ivif` | 红外–可见光论文最常报告的 8 个指标：EN（`Q_En`）、SD（`Q_CI`）、SF、MI（`Q_MI`）、SCD（`Q_SCD`）、VIF（`VIF_Sheikh`）、Qabf（`Q_ABF`）、SSIM（`MS_SSIM`） |
| `image` | 只依赖融合图像：`Q_MLI`、`Q_En`、`Q_CI`、`SF`、`AG_VIFB`、`Q_EI` |
| `information` | 信息论：`Q_MI`、`Q_NMI`、`Q_NCIE`、`Q_TE_Nava`、`Q_FMI`、`Q_CE` |
| `feature` | 边缘 / 特征：`Q_ABF`、`Q_M`、`Q_SF`、`Q_FSIM` |
| `structural` | 结构相似性：`MS_SSIM`、`Q_UIQI` |
| `similarity` | 与源图像的相似性：`Q_PSNR`、`Q_CC`、`Q_SCD`、`VIF_Sheikh`、`Q_VIFF` |

共 35 个指标。每个指标的原始论文、取值范围和示例见 [README](README.md#-metrics) 与 [API 文档](docs/API.md)。

## 常见问题：为什么我的数值对不上？

1. **尺度不一致。** 融合图在 [0,1]、源图像在 [0,255]（或反过来）时，Qabf 会接近 0：同一组图像在尺度一致时
   Qabf = 0.31，把融合图改成 [0,1] 后只有 0.0003。请把三张图统一到 0–255。
2. **灰度转换不同。** MATLAB `rgb2gray` 与 OpenCV 的系数不同，伪彩色 PET/SPECT 图像差异尤其明显。
   `fm.read_gray` 复现 MATLAB，`fm.read_gray_cv` 使用 OpenCV。
3. **MATLAB 特有细节。** 边界处理、四舍五入规则、`uint8` 运算饱和等，本库均按原版实现。

## 验证

```bash
python verify/selftest.py      # 输出 ALL OK
```

参考值由原版 `.m` 文件在 GNU Octave 中生成，测试用例包括覆盖各种边界情况的合成图像和公开的 Harvard 全脑图谱
PET/SPECT–MRI 图像对。详见 [README › Verification](README.md#-verification)。

## 路线图与贡献

还有 9 个指标在[路线图](https://github.com/VuVietAnh2309/fusion-metrics/issues?q=is%3Aopen+label%3Ametric-port)上。
欢迎提 issue 或 PR（可以用中文）；如果发现某个值与 MATLAB 不一致，请附上图像和 MATLAB 的结果。贡献指南见
[CONTRIBUTING.md](CONTRIBUTING.md)。

## 引用

如果本库对你的研究有帮助，请引用（DOI：[10.5281/zenodo.24536175](https://doi.org/10.5281/zenodo.24536175)），
BibTeX 见 [README › Citation](README.md#-citation)。也欢迎在 GitHub 上点个 ⭐。

## 许可证

MIT。部分由第三方工具箱移植的指标仍遵循原作者的条款（仅限教育与科研用途），详见 [NOTICE](NOTICE)。
