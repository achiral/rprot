# 分析・開発用のコード
#
# install_packages.R でインストールしたパッケージを普通にロードして使用する。
# 例:
#     library(palmerpenguins)
#     head(penguins)

setwd("/home/rstudio/project")
# パッケージはhome/rstudio/project/dev/packages.Rにて管理すると良い
# 新しいパッケージをライブラリにインストールするときは `renv::install()` 関数を使用
# renv::install("tidyverse")


# ライブラリの状態を記録するには `renv::snapshot()` 関数を使用
# ライブラリにインストールされたパッケージとそのバージョン情報を [`renv.lock`](./renv.lock) ファイルに記録
# GitHub にプッシュし、チームメンバー間で共有
# renv::snapshot()


# `renv.lock` ファイルのライブラリ状態を復元するときは`renv::restore()` 関数を使用
# renv::snapshot()