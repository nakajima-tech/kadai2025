# 2025年新入生課題　　　　　中嶋滉太

# 課題A

①ssh
②公開鍵
③authorized_keys
④755
⑤755
⑥644
⑦600
⑧秘密鍵
⑨600

# 課題B

1.CA:α炭素、CB:β炭素、CD:δ炭素、CG:γ炭素、C:カルボニル炭素、O:カルボニル酸素、N:主鎖アミド窒素
2.鎖識別子であり、A鎖に含まれることを示す
3.同一の原子に複数の存在確率があることを示している。この時、55~60文字目の値は0.50になっている。
4.温度(β)因子であり、原子の熱振動の大きさを表す。値が高いと原子の結晶内での動きは大きく、低いと小さいことを意味する。値が高いほど構造が不安定であり、低いほど安定となる。
5.結晶化剤：SO4,GOL,HOH
　リガンド：GIX

# 課題C

1.

![課題C-1](figures/kadaiC-1_1alk_zanki_stick.png)
2.ASP-51,ASP-101,SER102,ALA-103,ASP-153,ALA-154,ARG-166,ASP-327,LYS-328,HIS-331,ASP-369,HIS-370,HIS-412
3.

![課題C-3](figures/kadaiC-3_1alk_1ew2_hikaku.png)

4.

![課題C-4](figures/kadaiC-4_shinsuisei_sosuisei.png)
　タンパク質構造中において、疎水性アミノ酸はタンパク質の内部に親水性アミノ酸は外部表面に概して存在している。

5.

![課題C-5](figures/kadaiC-5.png)
保存度9となっている残基はタンパク質の内部や活性部位に分布している。

# 課題D

空欄のコード

```python:kadaiD
    Z = np.array([activation(np.array([x1, x2])) for x1, x2 in zip(xx1.ravel(), xx2.ravel())])
    Z = Z.reshape(xx1.shape)
```

学習回数30回、学習率0.05時のprint(errors_list),print(w)の結果
[19, 17, 13, 12, 11, 10, 8, 6, 5, 3, 4, 3, 4, 3, 2, 3, 2, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
[-12.8    1.62   1.59]

線形分離後の画像
![課題D](figures/kadaiD-5_gakusyugo_bunri.png)

# 課題E

1.19.8%
2.pythonファイル

# 課題F

github
ユーザ名：nakajima-tech   レポジトリ：kadaie1 で公開
URL:<https://github.com/nakajima-tech/kadaie1.git>

# 課題G

1.

```python:kadaiG-1
import numpy as np

def needleman_wunsch(seq1, seq2, match=5, mismatch=-2, gap=-6):
    nx = len(seq1)
    ny = len(seq2)
    F = np.zeros([nx + 1, ny + 1], dtype=int)

    # 1行目の値を埋める。
    for i in range(1, nx + 1):
        F[i, 0] = F[i - 1, 0] + gap

    # 1列目の値を埋める。
    for j in range(1, ny + 1):
        F[0, j] = F[0, j - 1] + gap

    # 最適なF[i,j]要素を計算してFの行列の要素を埋めていく
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            # ここにアルゴリズムを書く
            if seq1[i - 1] == seq2[j - 1]:
                score = match
            else:
                score = mismatch

            diag = F[i - 1, j - 1] + score
            up = F[i - 1, j] + gap
            left = F[i, j - 1] + gap
            F[i, j] = max(diag, up, left)

    # 擬似コードの部分。最適化された配列2つが返される。
    optseq1 = ""
    optseq2 = ""
    i = nx
    j = ny

    while i > 0 or j > 0:
        if i > 0 and F[i, j] == F[i - 1, j] + gap:
            optseq1 = seq1[i - 1] + optseq1
            optseq2 = "-" + optseq2
            i -= 1
        elif j > 0 and F[i, j] == F[i, j - 1] + gap:
            optseq1 = "-" + optseq1
            optseq2 = seq2[j - 1] + optseq2
            j -= 1
        else:
            optseq1 = seq1[i - 1] + optseq1
            optseq2 = seq2[j - 1] + optseq2
            i -= 1
            j -= 1
    # 最後にtuple型で最適化された2つの配列をreturnする。
    return optseq1, optseq2

seq1 = "CTAAGGGATTCCGGTAATTAGACAG"
seq2 = "ATAGACCATATGTCAGTGACTGTGTAA"

output1, output2 = needleman_wunsch(seq1, seq2)
print(output1 + "\n" + output2)
```

2.

```python:kadaiG-2
import numpy as np

dna_similarity_matrix = {
    "A": {"A": 10, "G": -1, "C": -3, "T": -4},
    "G": {"A": -1, "G": 7, "C": -5, "T": -3},
    "C": {"A": -3, "G": -5, "C": 9, "T": 0},
    "T": {"A": -4, "G": -3, "C": 0, "T": 8},
}

def s(A, B, match=None, mismatch=None):
    if (match is not None) or (mismatch is not None):
        if match is None:
            match = 5
        if mismatch is None:
            mismatch = -2
        return match if A == B else mismatch
    else:
        return dna_similarity_matrix[A][B]

def needleman_wunsch(seq1, seq2, match=None, mismatch=None, gap=-6):
    nx = len(seq1)
    ny = len(seq2)
    F = np.zeros([nx + 1, ny + 1], dtype=int)

    # 1行目の値を埋める。
    for i in range(1, nx + 1):
        F[i, 0] = F[i - 1, 0] + gap

    # 1列目の値を埋める。
    for j in range(1, ny + 1):
        F[0, j] = F[0, j - 1] + gap

    # 最適なF[i,j]要素を計算してFの行列の要素を埋めていく
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            # ここにアルゴリズムを書く
            score = s(seq1[i - 1], seq2[j - 1], match, mismatch)
            diag = F[i - 1, j - 1] + score
            up = F[i - 1, j] + gap
            left = F[i, j - 1] + gap
            F[i, j] = max(diag, up, left)

    # 擬似コードの部分。最適化された配列2つが返される。
    optseq1 = ""
    optseq2 = ""
    i = nx
    j = ny

    while i > 0 or j > 0:
        if i > 0 and F[i, j] == F[i - 1, j] + gap:
            optseq1 = seq1[i - 1] + optseq1
            optseq2 = "-" + optseq2
            i -= 1
        elif j > 0 and F[i, j] == F[i, j - 1] + gap:
            optseq1 = "-" + optseq1
            optseq2 = seq2[j - 1] + optseq2
            j -= 1
        else:
            optseq1 = seq1[i - 1] + optseq1
            optseq2 = seq2[j - 1] + optseq2
            i -= 1
            j -= 1
    # 最後にtuple型で最適化された2つの配列をreturnする。
    return optseq1, optseq2

seq1 = "CTTCACC"
seq2 = "GTTTCACT"

output1, output2 = needleman_wunsch(seq1, seq2)
print(output1 + "\n" + output2)
```

3.

```python:kadaiG-3
import numpy as np

aminodict = dict(
    A=0, R=1, N=2, D=3, C=4, Q=5, E=6, G=7, H=8, I=9,
    L=10, K=11, M=12, F=13, P=14, S=15, T=16, W=17, Y=18, V=19,
)

BLOSUM62_MATRIX = [
    [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
    [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
    [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
    [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
    [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
    [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
    [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
    [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
    [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
    [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
    [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
    [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
    [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
    [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
    [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
    [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
    [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
]

def score(a, b):
    return BLOSUM62_MATRIX[aminodict[a]][aminodict[b]]

def needleman_wunsch(seq1, seq2, gap=-6):
    nx, ny = len(seq1), len(seq2)
    F = np.zeros((nx + 1, ny + 1), dtype=int)

    for i in range(1, nx + 1):
        F[i, 0] = F[i - 1, 0] + gap
    for j in range(1, ny + 1):
        F[0, j] = F[0, j - 1] + gap


    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            match_score = score(seq1[i - 1], seq2[j - 1])
            F[i, j] = max(
                F[i - 1, j - 1] + match_score,
                F[i - 1, j] + gap,
                F[i, j - 1] + gap
            )


    optseq1, optseq2 = "", ""
    i, j = nx, ny
    total_score = 0
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i, j] == F[i - 1, j - 1] + score(seq1[i - 1], seq2[j - 1]):
            optseq1 = seq1[i - 1] + optseq1
            optseq2 = seq2[j - 1] + optseq2
            total_score += score(seq1[i - 1], seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and F[i, j] == F[i - 1, j] + gap:
            optseq1 = seq1[i - 1] + optseq1
            optseq2 = "-" + optseq2
            total_score += gap
            i -= 1
        else:
            optseq1 = "-" + optseq1
            optseq2 = seq2[j - 1] + optseq2
            total_score += gap
            j -= 1

    return optseq1, optseq2, total_score

if __name__ == "__main__":
    seq1 = input("Seq1: ").strip().upper()
    seq2 = input("Seq2: ").strip().upper()

    optseq1, optseq2, score_val = needleman_wunsch(seq1, seq2)

    print(optseq1)
    print(optseq2)
    print(f"\nScore: {score_val}")
```

# 課題H

1.RMSD =    0.500 (3157 to 3157 atoms)

重ね合わせの図
![課題H-1](figures/kadaiH-1_alphafold0001.png)

2.

A9FZ87,G5EKN0,B5GW45におけるそれぞれのOPPと1位炭素または3位炭素の距離、または環化が起きるはずの1位-10/11位炭素距離間, 1-6/7位炭素距離間の測定結果は下記のようになった。
OPPと1位または3位の距離はOPPを構成する全ての原子に対して行った。
A9FZ87
Distance C1 - PA: 2.65 Å
Distance C1 - PB: 4.42 Å
Distance C1 - O1A: 3.85 Å
Distance C1 - O2A: 3.21 Å
Distance C1 - O3A: 2.86 Å
Distance C1 - O1B: 5.23 Å
Distance C1 - O2B: 4.80 Å
Distance C1 - O3B: 5.06 Å　　　最大値：5.23 Å　最小値：2.65 Å

Distance C3 - PA: 4.41 Å
Distance C3 - PB: 5.45 Å
Distance C3 - O1A: 5.12 Å
Distance C3 - O2A: 5.49 Å
Distance C3 - O3A: 4.12 Å
Distance C3 - O1B: 5.84 Å
Distance C3 - O2B: 5.44 Å
Distance C3 - O3B: 6.59 Å　　　最大値：6.59 Å　最小値：4.12 Å

Distance C1 - C10: 6.99 Å
Distance C1 - C11: 5.59 Å
Distance C1 - C6: 4.87 Å
Distance C1 - C7: 4.55 Å

G5EKN0
Distance C1 - PA: 2.62 Å
Distance C1 - PB: 4.58 Å
Distance C1 - O1A: 3.84 Å
Distance C1 - O2A: 3.06 Å
Distance C1 - O3A: 3.02 Å
Distance C1 - O1B: 5.09 Å
Distance C1 - O2B: 5.16 Å
Distance C1 - O3B: 5.33 Å     最大値：5.33 Å　最小値：2.62 Å

Distance C3 - PA: 4.10 Å
Distance C3 - PB: 4.46 Å
Distance C3 - O1A: 5.35 Å
Distance C3 - O2A: 4.63 Å
Distance C3 - O3A: 3.33 Å
Distance C3 - O1B: 4.14 Å
Distance C3 - O2B: 5.26 Å
Distance C3 - O3B: 5.49 Å     最大値：5.49 Å　最小値：3.33 Å

Distance C1 - C10: 4.73 Å
Distance C1 - C11: 7.13 Å
Distance C1 - C6: 4.59 Å
Distance C1 - C7: 5.36 Å

B5GW45
Distance C1 - PA: 2.62 Å
Distance C1 - PB: 4.38 Å
Distance C1 - O1A: 3.78 Å
Distance C1 - O2A: 3.25 Å
Distance C1 - O3A: 2.83 Å
Distance C1 - O1B: 5.24 Å
Distance C1 - O2B: 4.80 Å
Distance C1 - O3B: 4.94 Å     最大値：5.24 Å　最小値：2.62 Å

Distance C3 - PA: 4.13 Å
Distance C3 - PB: 5.48 Å
Distance C3 - O1A: 4.62 Å
Distance C3 - O2A: 5.30 Å
Distance C3 - O3A: 4.06 Å
Distance C3 - O1B: 5.90 Å
Distance C3 - O2B: 5.63 Å
Distance C3 - O3B: 6.54 Å     最大値：6.54 Å　最小値：4.13 Å

Distance C1 - C10: 6.00 Å
Distance C1 - C11: 6.07 Å
Distance C1 - C6: 4.70 Å
Distance C1 - C7: 4.68 Å

# 課題I

1.チュートリアルを実行した
2.
![課題I-2](figures/kadaiI-2.png)
3.(1)剛体球(2)ばね(3)角度(4)二面角(5)座標(6)トポロジーファイル(7)系のエネルギー最小化(8)系のエネルギー平衡化(9)系のProduction Run(10)ミッシング領域(11)ジスルフィド(12)リガンド(13)antechamber

# 課題J

1.kadaiJ1-1のファイルをkadaiJ1-2ファイルと同じフォルダで実行することで動く
2.kadaiJ2-1のファイルをkadaiJ2-2ファイルと同じフォルダで実行することで動く

# 課題K

1.

フォルダ内のwordファイルに清書版をまとめました
**Cartesian Coordinates and Energies**

**IM5**

**B3LYP/6-31+G(d,p):**

**Sum of electronic and thermal Free Energies = -586.070532**

**Imaginary Frequencies : none**

C 0.865968 -0.907178 -1.029866

C 2.000013 -0.502431 0.029118

C 1.497891 0.513103 1.087495

H 0.997139 -0.248988 -1.899487

H 1.041009 -1.931188 -1.369672

H 2.371625 0.775167 1.695008

H 0.821065 -0.016565 1.777059

C -1.489526 -1.869214 -0.493331

C -0.532912 -0.73457 -0.613336

C -2.571931 0.150834 0.290303

C -2.881042 -1.216481 -0.333548

H -1.394249 -2.521697 -1.371989

H -2.426183 0.083220 1.373775

H -3.353112 0.888519 0.109531

H -3.555782 -1.825834 0.270938

H -3.347755 -1.085426 -1.316696

H -1.231889 0.545325 -0.393208

C 0.673030 2.332796 -0.740552

H 1.050984 1.751964 -1.572202

H 0.764812 3.402887 -0.900403

C -1.102805 -2.750884 0.742886

C -1.864243 -3.529875 0.827441

H -0.133177 -3.231959 0.610692

H -1.095973 -2.174676 1.671830

C 2.439503 -1.779500 0.775309

H 2.816780 -2.537544 0.081340

H 3.252074 -1.533994 1.466616

H 1.632759 -2.218995 1.367411

C 3.216688 0.034283 -0.750251

H 4.047923 0.202503 -0.058067

H 3.554287 -0.691064 -1.498287

C 3.014375 0.980513 -1.255740

C -1.531541 2.935463 0.462958

H -0.988267 3.857331 0.686429

H -2.260712 3.171685 -0.320552

H -2.077349 2.650382 1.367694

H -1.516529 0.656615 -1.480893

C -0.552639 1.857147 0.013657 C 0.817385 1.813869 0.671944

H 0.951338 2.591680 1.420758

**TS_5-6**

**B3LYP/6-31+G(d,p):**

**Sum of electronic and thermal Free Energies = -586.067179**

**Imaginary Frequencies : 1 (-601.77 cm^-1^)**

C 0.863719 -0.938850 -0.996132

C 2.002233 -0.555412 0.028129

C 1.568248 0.521418 1.060881

H 1.021911 -0.385317 -1.927675

H 0.966884 -1.992004 -1.272475

H 2.470120 0.791468 1.621608

H 0.906319 0.049170 1.801552

C -1.575295 -1.832876 -0.489233

C -0.560763 -0.707462 -0.571332

C -2.601122 0.257863 0.248796

C -2.936666 -1.108209 -0.374109

H -1.510133 -2.462415 -1.383907

H -2.587812 0.222621 1.345307

H -3.296156 1.047967 -0.038269

H -3.660456 -1.675693 0.213779

H -3.368274 -0.963062 -1.371541

C -1.168954 0.537386 -0.214074

C 0.702763 2.265470 -0.802913

H 1.024062 1.637436 -1.623483

H 0.776171 3.325580 -1.025459

C -1.258483 -2.724509 0.737579

H -2.023212 -3.503004 0.800571

H -0.288940 -3.217940 0.645310

H -1.275110 -2.157297 1.673510

C 2.402336 -1.814094 0.827820

H 2.761280 -2.606860 0.163142

H 3.213366 -1.576439 1.523499

H 1.571746 -2.210759 1.419496

C 3.236252 -0.090515 -0.774164

H 4.085743 0.063452 -0.101510

H 3.532714 -0.847950 -1.508023

H 3.063572 0.848150 -1.307298

C -1.443586 2.983651 0.467292

H -0.859914 3.887428 0.658706

H -2.179019 3.231907 -0.305746

H -1.978621 2.741129 1.390806

H -1.204359 0.264400 -1.425500

C -0.514440 1.859681 0.025637

C 0.910923 1.822414 0.611564

H 1.075131 2.632296 1.319170

**IM6**

**B3LYP/6-31+G(d,p):**

**Sum of electronic and thermal Free Energies = -586.087938**

**Imaginary Frequencies : none**

C 0.994770 -0.824020 -1.059527

C 2.041753 -0.388853 0.025524

C 1.506368 0.635822 1.057719

H 1.232981 -0.286500 -1.982428

H 1.184439 -1.874022 -1.307811

H 2.372624 1.030583 1.605441

H 0.905545 0.124723 1.820124

C -1.365832 -1.976624 -0.424589

C -0.544714 -0.713165 -0.862286

C -2.612546 0.014981 0.189095

C -2.798730 -1.407804 -0.365940

H -1.281600 -2.726774 -1.217073

H -2.629429 0.019064 1.292482

H -3.380516 0.736401 -0.110441

H -3.472584 -2.007796 0.249179

H -3.231148 -1.356596 -1.371690

C -1.236203 0.447609 -0.213100

C 0.472767 2.288626 -0.763949

H 0.817192 1.695930 -1.600079

H 0.393529 3.352017 -0.970287

C -0.933495 -2.624894 0.895259

H -1.627534 -3.431712 1.155790

H 0.060693 -3.068785 0.819137

H -0.925124 -1.915151 1.732919

C 2.551438 -1.602493 0.833370

H 2.980886 -2.359062 0.168749

H 3.338144 -1.294471 1.530426

H 1.761314 -2.073425 1.422225

C 3.261233 0.206219 -0.718790

H 4.070108 0.433702 -0.017101

H 3.652244 -0.509174 -1.449922

H 3.019380 1.128573 -1.256857

C -1.765923 2.819660 0.508430

H -1.242245 3.758840 0.701750

H -2.550148 3.030636 -0.226332

H -2.245445 2.511763 1.442578

H -0.939460 -0.543435 -1.890231

C -0.801751 1.766066 -0.012386

C 0.739656 1.857104 0.595499

H 0.683398 2.635277 1.353482

2.

活性化エネルギー：2.10 kcal/mol

3.

虚振動が大きい場合は小さい場合と比べて遷移状態の構造が不安定であるということを表し、虚振動が大きい場合ほど反応物と遷移状態の間のエネルギーダイヤグラムの差は小さくなる。

4.

PM7とB3LYPで計算した構造間には差が見られた。IM5とTS状態においては3員環と5員環に結合する炭素原子の飛び出す方向が異なっている。IM6においては違いは見られなかった。

5.

今回の遷移状態構造は出発物の構造に近い。
これは、遷移状態と反応物、生成物とのエネルギー差が近いほど構造状態も近くなるというハモンドの仮説より説明できる。今回のSum of electronic and thermal Free EnergiesはIM5が-586.070532、IM5-6が-586.067179、IM6は-586.087938であり、IM5とIM5-6間のエネルギー差は0.003353であり、IM6とIM5-6間は0.020759であった。よって、今回は出発物と遷移状態とのエネルギー差の方が小さいため、構造も類似したと考えられる。

# 課題L

完成させたガウス過程回帰実装のコードを下記に示す。

```python:kadaiL

import matplotlib.pyplot as plt
import numpy as np
from typing import Tuple

rng = np.random.default_rng(3)


def func(x):
    return np.sin(2 * np.pi * x) + 0.3 * x


def create_toy_data(func, low=0, high=1.0, n=10, std=0.1):
    x_train = rng.uniform(low, high, n)
    t_train = func(x_train) + rng.normal(scale=std, size=n)
    return x_train, t_train


x_train = np.array(
    [
        0.08564917,
        0.23681051,
        0.80127447,
        0.58216204,
        0.09412864,
        0.43312694,
        0.4790513,
        0.15973891,
        0.73457715,
        0.11367202,
    ]
)
t_train = np.array(
    [
        0.56082141,
        1.03234815,
        -0.73629794,
        -0.38576905,
        0.48027033,
        0.49877898,
        0.32315477,
        0.86751411,
        -0.67915941,
        0.66915142,
    ]
)

fig, ax = plt.subplots(dpi=300, figsize=(5, 5))
ax.scatter(x_train, t_train, color="blue", label="observation")
ax.plot(
    np.linspace(0, 1, 100),
    func(np.linspace(0, 1, 100)),
    color="blue",
    label="ground truth",
)
ax.set_xlabel("$x$")
ax.set_ylabel("$t$")
ax.legend()
plt.show()


class MyKernel:
    def __init__(self, thetas):
        assert np.shape(thetas) == (4,)
        self.thetas = thetas

    def get_thetas(self):
        return np.copy(self.thetas)

    def get_sqdist(self, x, y):
        x = np.atleast_2d(x)
        y = np.atleast_2d(y)
        sqdist = (
            np.sum(x**2, axis=1)[:, None] + np.sum(y**2, axis=1)[None, :] - 2 * x @ y.T
        )
        return sqdist

    def __call__(self, x, y):
        x = x[:, None]
        y = y[:, None]
        sqdist = self.get_sqdist(x, y)
        theta0, theta1, theta2, theta3 = self.thetas
        K = theta0 * np.exp(-theta1 / 2 * sqdist) + theta2 + theta3 * x @ y.T
        return K

    def derivatives(self, x, y):
        x = x[:, None]
        y = y[:, None]
        sqdist = self.get_sqdist(x, y)
        theta0, theta1, theta2, theta3 = self.thetas
        delta_0 = np.exp(-theta1 / 2 * sqdist)
        delta_1 = -0.5 * theta0 * sqdist * np.exp(-theta1 / 2 * sqdist)
        delta_2 = np.ones_like(sqdist)
        delta_3 = x @ y.T
        return (delta_0, delta_1, delta_2, delta_3)

    def update_parameters(self, updates):
        assert np.shape(updates) == (4,)
        self.thetas += updates


class MyGPR:
    def __init__(self, kernel, beta=1.0):
        self.kernel = kernel
        self.beta = beta

    def fit_kernel(self, x, t, learning_rate=0.1, iter_max=10000):
        self.x = x
        self.t = t
        for i in range(iter_max):
            thetas = self.kernel.get_thetas()
            K = self.kernel(x, x)
            self.covariance = K + np.identity(x.shape[0]) / self.beta
            self.precision = np.linalg.inv(self.covariance)
            gradients = self.kernel.derivatives(x, x)
            updates = np.array(
                [
                    0.5
                    * np.trace(
                        (
                            self.precision @ (t[:, None] @ t[None, :]) @ self.precision
                            - self.precision
                        )
                        @ grad
                    )
                    for grad in gradients
                ]
            )
            self.kernel.update_parameters(learning_rate * updates)
            if np.allclose(thetas, self.kernel.get_thetas()):
                break
        else:
            print("parameters may not have converged")

    def predict(self, x):
        x = x[:, None]
        kappa = self.kernel(self.x, x.ravel())
        mean = kappa.T @ self.precision @ self.t
        k_star = self.kernel(x.ravel(), x.ravel())
        std = np.sqrt(np.maximum(0, np.diag(k_star - kappa.T @ self.precision @ kappa)))
        return mean, std


kernel = MyKernel(thetas=np.array([0.5, 4.0, 0.0, 0.0]))
K = kernel(x_train, x_train)
print(f"K = {K}")

regression = MyGPR(kernel=kernel, beta=100)
regression.fit_kernel(x_train, t_train, learning_rate=0.1, iter_max=5000)

x_test = np.linspace(0, 1, 100)
y, y_std = regression.predict(x_test)

fig, ax = plt.subplots(dpi=300, figsize=(5, 5))
ax.scatter(x_train, t_train, alpha=0.5, color="blue", label="observation")
ax.plot(x_test, func(x_test), color="blue", label="ground truth")
ax.plot(x_test, y, color="red", label="predict_mean")
ax.fill_between(
    x_test, y - y_std, y + y_std, color="pink", alpha=0.5, label="predict_std"
)
ax.legend(loc="lower left")
ax.set_xlabel("$x$")
ax.set_ylabel("$t$")
ax.legend()
plt.show()
```

# 課題M

![課題M](figures/kadaiM.png)
