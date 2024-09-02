## WARNING

**This package is experimental and not supported. To fit a non-parametric model for causal inference, use [bartCause](https://github.com/vdorie/bartCause) instead.**

## Prequisites

[Rtools](https://cran.r-project.org/bin/windows/Rtools/) for windows. [gfortran](https://github.com/fxcoudert/gfortran-for-macOS/releases) for Mac OS.

## Notes

Contains R code to fit non-parmetric causal inference models. Of particular interest are the calculations in the examples folder.

To run the examples, first install the package by running, from within R

    install.packages("path/to/repository", type = "source", repos = NULL)

or from the command line

    R CMD INSTALL path/to/repository

The file `examples/toy/generatePlots.R` contains code to produce the toy data example. The simulation data can be simply loaded by running `examples/ihdp_sim/loadData.R`, or the simulations produced by running either `examples/ihdp_sim/runLocally.R` or `examples/ihdp_sim/queueJobs.R` (if on a cluster using TORQUE). Change your working directory before doing any of the above.

<!-- ###############################################################################################################
## 警告

**このパッケージは実験的なものであり、サポートされていません。因果推論のためにノンパラメトリックモデルをフィットするには、代わりに [bartCause](https://github.com/vdorie/bartCause) を使用してください。

## 前提条件

Windows用の[Rtools](https://cran.r-project.org/bin/windows/Rtools/). [gfortran](https://github.com/fxcoudert/gfortran-for-macOS/releases) Mac OS用。

## ノート

ノンパラメトリック因果推論モデルのRコード。特に興味深いのは、examplesフォルダにある計算です。

サンプルを実行するには、まずRからパッケージをインストールしてください。

    インストール.packages("path/to/repository", type = "source", repos = NULL)

またはコマンドラインから

    R CMD INSTALL path/to/repository

examples/toy/generatePlots.R`ファイルには、トイデータの例を作成するコードが含まれています。シミュレーションデータは、`examples/ihdp_sim/loadData.R`を実行してロードするか、`examples/ihdp_sim/runLocally.R`または`examples/ihdp_sim/queueJobs.R`（TORQUEを使用するクラスタの場合）を実行して生成します。上記のいずれかを行う前に、作業ディレクトリを変更してください。 -->
