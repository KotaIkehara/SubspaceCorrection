# Automated Subspace Correction法を前処理とするCGソルバ
以下では，不完全コレスキー分解（IC分解）を前処理とする共役勾配法（CG法）に基づく線型ソルバを***ICCG***，Subspace Correction法[[1]](https://pdfs.semanticscholar.org/edf9/7016849f3def445b07d33a6c30649a5036c6.pdf)を前処理とするICCG法を***SCICCG***と表記する．

## ファイルの説明
- main.iccg.c：逐次版ICCG
-  main.iccg-OpenMP.c：並列版ICCG
-  main.sciccg.c：逐次版SCICCG
-  main.sciccg-OpenMP.c：並列版SCICCG

並列版と記載のあるソルバでは，[OpenMP](https://www.openmp.org/)を用いたマルチスレッド並列化を行っている．
