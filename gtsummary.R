# はじめに

# Rで統計解析を行う際、Rで結果がわかっても、患者の背景表(table1)や結果を論文やスライドに転記する際手入力だと時間がかかり再現性の問題もでてきます。 それまではtable1に関しては日本人の先生が作ったtableoneパッケージ（https://cran.r-project.org/web/packages/tableone/index.html）が主流でした。 しかし、2020年あたりにできたgtsummaryは、table1だけでなく、回帰分析やロジスティック回帰分析などの結果をわかりやすく出力できることで、利用者が一気に増えました。

# 今回はmJHONSNOW内でも紹介されているgtsummaryが使えるようになることを目標にします。

## 困ったらここを見よう

# 公式：<https://www.danieldsjoberg.com/gtsummary/index.html>

# CursorやChatGPTにこのサイトを引用して解説してもらうだけでも大抵のことはわかると思います。

# gtsummaryができること

# 1.  背景表（table1）を作る　→　tbl_summary()
# 2.  多変量解析の結果を出力する　→　tbl_regression()
# 3.  2×2の表を出力する(今回は割愛)　→　tbl_cross()

# 第1章 `tbl_summary()` が使えるようになる

## 0. パッケージとデータの準備

# まずは今回使うパッケージを読み込んでおきましょう。

# 今回はpacmanパッケージを使ってパッケージを読み込みます。pacmanパッケージは`install.packages()`と`library()`を同時に読み込むことができるパッケージです。

#pacmanパッケージが入っていない方はインストールする

if (!require("pacman")) install.packages("pacman")

#残りのパッケージを読み込む
pacman::p_load(
    gtsummary, #今回の本体
    tidyverse, #データの前処理
    broom, #統計の結果をデータフレームに変換する
    broom.helpers, #統計の結果をデータフレームに変換する
    gt, #必須ではないが、細かいカスタマイズが可能
    flextable, #必須ではないが、細かいカスタマイズが可能
    survival, #pbcデータセットが入っている
    causaldata, #重回帰分析で使うデータセットが入っている
    MatchIt, #傾向スコアマッチング
    renv, #再現性のためのパッケージ
    webshot2, #gtsummaryの結果を画像として保存するためのパッケージ
    car, #VIFのためのパッケージ
    cardx, #VIFのためのパッケージ
    smd, #SMDの計算
    ggstats, #回帰係数のグラフで必要,
    knitr,
    rmarkdown,
    languageserver #cursor使う人向け
)

# 今回使うパッケージの私のバージョン。今回の勉強会向けのためだけなので、以下は不要
si <- sessionInfo()
# 基本パッケージ以外でロード済みのパッケージ一覧
loaded_pkgs <- c(
    "gtsummary",
    "tidyverse",
    "broom",
    "broom.helpers",
    "gt",
    "flextable",
    "survival",
    "causaldata",
    "MatchIt",
    "renv",
    "webshot2",
    "car",
    "cardx",
    "smd",
    "ggstats"
)
# バージョンを抽出
versions <- sapply(loaded_pkgs, function(pk) {
    as.character(packageVersion(pk))
})
# データフレーム化
df <- data.frame(
    Package = loaded_pkgs,
    Version = versions,
    row.names = NULL,
    stringsAsFactors = FALSE
)
print(df)

# もし今回のコードの実行中以下のようなエラーが出たら**パッケージのバージョンが古い**ので、パッケージを更新する必要があります。

# たとえばこの図ではbroom.helpersパッケージのエラーなので、`install.packages("broom.helpers")`でパッケージを更新します。

# ここではsurvivalパッケージの `pbc`（原発性胆汁性肝硬変）データを使って基本的な使い方を確認していきます。

p_load(survival) # 本当は実行済みだが勉強のため記載
data(pbc) # survivalパッケージの中にあるpbcというdata.frameを読み込んでという意味

# | 変数名 | クラス | 説明 |
# |------------------------|------------------------|------------------------|
# | `id` | integer | 症例番号 |
# | `time` | integer | 登録から死亡、移植、または1986年7月の研究分析までの日数 |
# | `status` | integer | エンドポイント時の状態, 0/1/2 は censored, transplant, dead |
# | `trt` | integer | 治療群（1=D-penicillmain, 2=placebo, NA=not randomised） |
# | `age` | numeric | 年齢（歳） |
# | `sex` | factor | 性別（m=男性, f=女性） |
# | `ascites` | integer | 腹水の有無 |
# | `hepato` | integer | 肝腫大の有無 |
# | `spiders` | integer | くも状血管腫の有無 |
# | `edema` | numeric | 浮腫スコア（0=なし, 0.5=未治療または治療成功, 1=利尿療法にもかかわらず浮腫） |
# | `bili` | numeric | 血清ビリルビン (mg/dL) |
# | `chol` | integer | 血清コレステロール (mg/dL) |
# | `albumin` | numeric | 血清アルブミン (g/dL) |
# | `copper` | integer | 尿中銅 (ug/day) |
# | `alk.phos` | numeric | アルカリホスファターゼ (U/liter) |
# | `ast` | numeric | AST（GOT） (U/ml) |
# | `trig` | integer | トリグリセリド (mg/dL) |
# | `platelet` | integer | 血小板数 |
# | `protime` | numeric | 標準化された血液凝固時間 |
# | `stage` | integer | 組織学的病期 (生検が必要) |

# 事前にデータの前処理をしておきます。

# trtのNAをfilterで除外
pbc <- pbc |> filter(!is.na(trt))

head(pbc)

## 0.1 表1を作る前に考えておくべきこと

# 表1（Table 1）を作成する前に、以下の3つの重要な点について検討する必要があります。

### 1. 変数の選択

# まず、どの変数を表1に含めるかを慎重に検討します。一般的に以下の変数は除外することが多いです：

# -   症例番号（ID）などの識別子

# -   主要アウトカム変数

# -   解析に使用しない変数

# 変数選択の際は、読者が研究対象集団の特徴を理解できるよう、臨床的に重要な変数を選択します。

### 2. 変数の尺度水準の決定

# 各変数の尺度水準（測定尺度）をどのように扱うかを決定します：

# -   連続尺度（例：年齢、血圧）

# -   順序尺度（例：重症度分類）

# -   カテゴリー変数（例：性別、治療群）

# 特に注意が必要なのは：

# -   2値データ（0/1）は通常カテゴリー変数として扱う

# -   5段階評価などの順序尺度は、連続変数として扱うか、カテゴリー変数として扱うか検討が必要

# この決定により、適用される統計手法や表示方法が大きく変わります

### 3. 要約統計量の選択

# 連続変数の要約方法を決定します：

# -   平均値±標準偏差

# -   中央値 \[四分位範囲\]

# 選択の基準：

# -   データの分布（正規分布しているか）

# -   外れ値の有無

# -   対象疾患での慣習

# -   先行研究での報告形式

# これらの決定は、データの特性を確認しながら行い、臨床的な意義と統計学的な適切性の両方を考慮する必要があります。`tableone`パッケージも`gtsummary`パッケージも、これらの設定に基づいて適切な表を作成することが重要です。

# 今回は以下の設定とします。

# -   使わない変数：id, status
# -   層別化する変数：trt
# -   数値(continuous)：time, age, bili, chol, albumin, copper, alk.phos, ast, trig, protime
# -   カテゴリ変数(categorical)：sex, ascites, hepato, spiders, edema, stage

# その上で数値は平均±標準偏差、で集計するものと中央値(四分位範囲)で集計するものを選択します。 - 平均値セット：time, age, albumin - 中央値セット：bili, chol, copper, alk.phos, ast, trig, protime

# 変数の設定を行うときはdput()が便利です。この結果をコピペすると変数名の選択が少し楽になるかもしれません。

# names()を使うと列名はわかる
names(pbc)
# この結果にdput()を使うとベクトル表記で返ってくる
names(pbc) |> dput()

# 変数設定
# varsと書いてあるのはvariablesの略で、列名の略称みたいなもの。別にvarsでなくてもいい。

# 除外する変数(必ずしも必要ではないが便利なこともある))
exclude_vars <- c("id", "status")
# 使う変数（この並びが行の順番になる）
include_vars <- c(
    "time",
    "age",
    "albumin",
    "bili",
    "chol",
    "copper",
    "alk.phos",
    "ast",
    "trig",
    "protime",
    "sex",
    "ascites",
    "hepato",
    "spiders",
    "edema",
    "stage"
)
# 数値型の変数
continuous_vars <- c("time", "age", "albumin")
median_vars <- c("bili", "chol", "copper", "alk.phos", "ast", "trig", "protime")
# カテゴリ変数
categorical_vars <- c("sex", "ascites", "hepato", "spiders", "edema", "stage")

# またもしp値ではなくSMDを表記する場合はfactor型に変換します。factor型にしないとtableoneパッケージと違うSMDを計算することを観測しました。そのためカテゴリー変数はとりあえずfactor型にしておくことをおすすめします。

# factor型にする列は今定義した`categorical_vars`です。ある列を作成・修正する時は`mutate`を使いますが、6列を1つずつ修正することになります。それでも構いませんが、今回はacrossを使いまとめて修正します。

# またtrtの列は層別化する列なのでfactor型にした上でラベルも作成します

# categorical_varsの列を全てfactor型に変換
pbc <-
    pbc |>
    mutate(
        across(all_of(categorical_vars), ~ factor(.x)),
        trt = factor(
            trt,
            levels = c(1, 2),
            labels = c("D-penicillmain", "placebo")
        ),
        trt = fct_rev(trt) # 後でロジスティック回帰を行う時に介入群を1とするためfactorのレベルを逆にする
    )

## 1. まずは触ってみよう

# まずは何も設定せず出力してみましょう。

# tbl_summary()は、デフォルトでは連続変数について中央値(第一四分位, 第三四分位)を表示します。

# 欠損を含むsummary
pbc |>
    tbl_summary()

# 次に、群ごとの要約表を作ってみましょう。`by =` 引数を使うことで、例えば治療群ごとの集計が可能になります。

# 治療群ごとの集計
pbc |>
    tbl_summary(by = trt)

# ここでid, status, status_factorは表示しないようにします。

# 表示しないようにするにはgtsummary内では`include =` 引数を使います。

# 表示する行の選択
pbc |>
    tbl_summary(
        by = trt,
        include = include_vars
    )

# デフォルトのgtsummaryは少し行の幅が広く、もし狭くするにはthemeを変えると変更できます。

# コンパクト表示にする
theme_gtsummary_compact()
pbc |>
    tbl_summary(
        by = trt,
        include = include_vars
    )

# 日本語のフォントを使う場合は`theme_gtsummary_language("ja")`を使います。

# 日本語のフォントを使う
theme_gtsummary_language("ja")
pbc |>
    tbl_summary(by = trt, include = include_vars)

# 雑誌に合わせた出力も可能です。

# 雑誌に合わせた出力
theme_gtsummary_journal("jama")
pbc |>
    tbl_summary(by = trt, include = include_vars) |>
    add_p()

# 一旦themeを戻す
reset_gtsummary_theme()

# p値を加えるには`add_p()`を使います。

# データの種類によって自動的に検定を選んでくれます。

# trtで群分けしp値を付ける
pbc |>
    tbl_summary(
        by = trt,
        include = include_vars
    ) |>
    add_p()

# また、gtsummaryのデフォルトは中央値（第1四分位, 第3四分位）です。 もし平均値（標準偏差）にする場合はthemeを変えると変更できます。

# 平均・標準偏差にする
theme_gtsummary_mean_sd()
pbc |>
    tbl_summary(
        by = trt,
        include = include_vars
    ) |>
    add_p()

# もしthemeをもとに戻す場合は`reset_gtsummary_theme`を使います。

# themeをリセットする
reset_gtsummary_theme()
pbc |>
    tbl_summary(
        by = trt,
        include = include_vars
    ) |>
    add_p()

## 2. いろいろな形で出力ができる

# 作成した表は、WordやPowerPoint、Excel形式など様々な形で出力できます。これは報告書やプレゼン資料を作る上で非常に重宝します。

library(gt)
library(flextable)

# themeをコンパクトにする
theme_gtsummary_compact()
summary_table <- pbc |>
    tbl_summary(
        by = trt,
        include = include_vars
    ) |>
    add_p() |>
    bold_labels()

# PowerPoint
summary_table |>
    as_flex_table() |>
    flextable::save_as_pptx(path = "summary_table.pptx")

# Word
summary_table |>
    as_flex_table() |>
    flextable::save_as_docx(path = "summary_table.docx")

# Excel形式（tibble）
as_tibble(summary_table) |>
    readr::write_excel_csv("sumamry_table.csv")

# 画像として保存
summary_table |> as_gt() |> gt::gtsave("summary_table.png")

# テキストとして出力
summary_table |> as_kable()


## 3. 細かいカスタムをしてみよう

# gtsummaryでは細かいカスタムをすることができます。

### 3.1 tbl_summaryのカスタマイズ

# tbl_summary()内で以下の調整を行うことができます。

# -   ラベルのカスタマイズ(label)
# -   集計方法のカスタマイズ(type)
# -   統計のカスタマイズ(statistic)
# -   有効数字のカスタマイズ(digits)

# どれも `選択した列名 ~ コマンド` の形で指定します。

### 3.1.1 ラベルのカスタマイズ

# コンパクト表示にする
theme_gtsummary_compact()

# ラベルをカスタマイズ1
pbc |>
    tbl_summary(
        by = trt,
        include = include_vars,
        label = age ~ "年齢 (歳)"
    )

# `list()`を使うと復数の設定が可能です。

# ラベルをカスタマイズ2
pbc |>
    tbl_summary(
        by = trt,
        include = include_vars,
        label = list(
            age ~ "年齢 (歳)",
            bili ~ "ビリルビン",
            albumin ~ "アルブミン"
        )
    )

### 3.1.2 集計方法のカスタマイズ

# 整数型の集計方法はここで選択することも可能です。

# -   連続変数：continuous
# -   カテゴリ変数：categorical
# -   二値変数：dichotomous

# ただ既にカテゴリーになっている変数はここでは設定できません。

# そこで**stageを整数型(integer)にしたら数値として認識するか？**を試してみます。

# タイプのカスタマイズ
pbc |>
    mutate(stage = as.integer(stage)) |>
    tbl_summary(
        by = trt,
        include = include_vars
    )

# しかし結果を見るとstageを整数型にしたのにカテゴリーとして集計されています。整数型はgtsummaryではカテゴリーとして集計されることが多いです。

# せ順序尺度として集計するには一旦整数型(integer)に変換した上で`type`を指定します。

# タイプのカスタマイズ2
pbc |>
    mutate(stage = as.integer(stage)) |>
    tbl_summary(
        by = trt,
        include = include_vars,
        type = list(stage ~ "continuous")
    )

# 他にも性別のような2値変数をもっと簡潔に表示したいという場合があります。

# そこで新たにmaleという列(男性だと1)を作ってみます。

# 2値変数をもっと簡潔に表示
pbc |>
    mutate(male = if_else(sex == "m", 1L, 0L)) |>
    tbl_summary(
        by = trt,
        include = c("male")
    )

# maleのような0/1で表現される列は二値変数：dichotomousとして認識され、1となった数と割合が表記されるようになります。

### 3.1.3 統計量のカスタマイズ

# 統計量のカスタマイズは`statistic`で行います。

# gtsummaryでは統計量を`{統計量}`の形で指定します。{}とその中に入っているもの以外はそのまま表記されます。

# 統計のカスタマイズ
pbc |>
    tbl_summary(
        by = trt,
        include = include_vars,
        statistic = list(
            continuous_vars ~ "{mean} ± {sd}",
            median_vars ~ "{median} [{IQR}]",
            categorical_vars ~ "{n} / {N} ({p}%)"
        )
    )

# `{統計量}`で使えるのは`{mean}`, `{sd}`, `{median}`, `{IQR}`, `{min}`, `{max}`, `{p◯}`, `{N}`などがあります。

# 応用編として、1つの列に2列の統計量など載せる方法もあります。

# その場合はtypeを`continnuous2`に設定することで解決できます。

pbc |>
    tbl_summary(
        by = trt,
        include = include_vars,
        type = c(continuous_vars, median_vars) ~ "continuous2",
        statistic = c(continuous_vars, median_vars) ~
            c(
                "{N_nonmiss}", # 欠損ではない数
                "{mean} ({sd})", # 平均(標準偏差)
                "{median} [{p25}, {p75}]"
            ), # 中央値[25%, 75%]
        missing = "no" # 欠損値の数は表示しない
    )

### 3.1.4 有効数字のカスタマイズ

# 有効数字を指定するにはdigitsを使います。

# 有効数字のカスタマイズ
pbc |>
    tbl_summary(
        by = trt,
        include = include_vars,
        digits = list(
            age ~ 5,
            albumin ~ 2,
            stage ~ 7
        )
    )

# 平均と標準偏差のように統計量が2つある場合、それぞれに有効数字をカスタマイズすることも可能です。

# stageの人数の有効桁数を0, ％を7桁にします。

# 有効数字のカスタマイズ2
pbc |>
    tbl_summary(
        by = trt,
        include = include_vars,
        digits = list(
            age ~ 5,
            albumin ~ 2,
            stage ~ c(0, 7)
        )
    )

### 3.1.15 欠損値の扱い

# 欠損値の指定はmissing = "ifany", "no", "always"の3択の選択肢があります。何も指定しなければifanyです。

pbc |>
    tbl_summary(
        by = trt,
        include = include_vars,
        missing = "always"
    )

# Unknownの文字を変えたい場合は missing_text = "◯◯"で変更できます

pbc |>
    tbl_summary(
        by = trt,
        include = include_vars,
        missing = "always",
        missing_text = "欠損値の数"
    )

### 3.2 別の情報の追加

# N_overallの追加
pbc |>
    tbl_summary(
        by = trt,
        include = include_vars
    ) |>
    add_n() |>
    add_overall() |>
    bold_labels()

### 3.3 検定の選択

# 検定はデータの種類により自動的に選択されます。特に連続変数に対しては中央値セットがデフォルトになります。

# もしデフォルトを平均値セットにするには`theme_gtsummary_mean_sd()`を使います。

# また中央値セットに戻したい場合は `reset_gtsummary_theme()` を使います。

# 平均値セットにする
theme_gtsummary_mean_sd()
pbc |>
    tbl_summary(by = trt, include = include_vars) |>
    add_p()
reset_gtsummary_theme()
theme_gtsummary_compact()

# 検定を1つ1つ決める場合は`add_p()`の中で`test =` 引数を使います。

# > -   "t.test" for a t-test,（t検定）
# >
# > -   "aov" for a one-way ANOVA test,（1元配置分散分析）
# >
# > -   "wilcox.test" for a Wilcoxon rank-sum test,（マン・ホットニーのU検定）
# >
# > -   "kruskal.test" for a Kruskal-Wallis rank-sum test,（クラスカル・ウォリス検定）
# >
# > -   "chisq.test" for a chi-squared test of independence,（χ二乗検定）
# >
# > -   "fisher.test" for a Fisher's exact test,（フィッシャーの正確確率検定）

# 最初に決めた設定を使い集計と検定を具体的に設定していきます。

# -   使わない変数：id, status
# -   層別化する変数：trt
# -   数値(continuous)：time, age, bili, chol, albumin, copper, alk.phos, ast, trig, protime
# -   カテゴリ変数(categorical)：sex, ascites, hepato, spiders, edema, stage

# その上で数値は平均(標準偏差)、で集計するものと中央値\[IQR\]で集計するものを選択します。

# -   平均値セット：time, age, albumin

# -   中央値セット：bili, chol, copper, alk.phos, ast, trig, protime

# -   カテゴリー変数：Fisherの正確確率検定

# 検定の選択
pbc |>
    tbl_summary(
        by = trt,
        include = include_vars,
        statistic = list(
            continuous_vars ~ "{mean}({sd})",
            median_vars ~ "{median}[{IQR}]"
        ),
        digits = list(
            continuous_vars ~ 1,
            median_vars ~ 1
        )
    ) |>
    add_p(
        test = list(
            continuous_vars ~ "t.test",
            median_vars ~ "wilcox.test",
            categorical_vars ~ "fisher.test"
        )
    )

### 3.4 SMDの出力

# 傾向スコアなどで表示するSMDは、最近は観察研究全般でもp値の代わりに表記されることが多いです。

# gtsummaryでは`add_difference()`を使うとSDMを出力できます。95%信頼区間も表示されます。

# 注意点としてはedemaのような数値型のカテゴリー変数はあらかじめfactor型にしておく必要があります。もししていない場合は数値としてのSDMが出力されます。

# SDMの出力
pbc_summary <- pbc |>
    tbl_summary(
        by = trt,
        include = include_vars,
        statistic = list(
            continuous_vars ~ "{mean}({sd})",
            median_vars ~ "{median}[{IQR}]"
        ),
        digits = list(
            continuous_vars ~ 1,
            median_vars ~ 1
        )
    )

pbc_summary |>
    add_difference(everything() ~ "smd") |>
    modify_column_hide("conf.low") # 入れなかったら95%信頼区間が表示される

### 3.5 ヘッダー、footnoteの調整

# ヘッダー、footnoteの調整は`modify_header()`、`modify_spanning_header()`、`modify_footnote()`を使います。

show_header_names(pbc_summary)

# ヘッダー、footnoteの調整
pbc_summary |>
    add_p() |>
    bold_labels() |>
    modify_header(
        label = "変数",
        all_stat_cols() ~ "**{level}<br>N = {n}({style_percent(p)}%)**"
    ) |>
    modify_spanning_header(all_stat_cols() ~ "マッチング前")

# 第2章 多変量解析の結果を出力しよう

## 1. 重回帰分析（nhefsデータ）

# 重回帰分析はアウトカムが連続変数なので、causaldataパッケージのnhefsデータを使います。

# nhefsデータの読み込み
nhefs <- causaldata::nhefs |>
    select(
        qsmk,
        sex,
        race,
        age,
        education,
        smokeintensity,
        smokeyrs,
        exercise,
        active,
        wt71,
        wt82_71
    )
head(nhefs)

# causaldataパッケージのnhefsデータセットの各変数の説明を作成
var_nhefs_desc <- tibble::tibble(
    変数 = c(
        "qsmk", # 禁煙指標（1=禁煙, 0=継続喫煙）
        "sex", # 性別（male=男性, female=女性）
        "race", # 人種（white=白人, black=黒人）
        "age", # ベースライン時の年齢（歳）
        "education", # 学歴（学校教育年数）
        "smokeintensity", # ベースライン時の1日あたり喫煙本数
        "smokeyrs", # ベースライン時の喫煙年数
        "exercise", # 運動レベル（カテゴリ変数）
        "active", # 身体活動の有無（カテゴリ変数）
        "wt71", # ベースライン時の体重（ポンド）
        "wt82_71" # 1971年から1982年の体重変化（ポンド）
    ),
    説明 = c(
        "禁煙指標（1=禁煙, 0=継続喫煙）",
        "性別（male=男性, female=女性）",
        "人種（white=白人, black=黒人）",
        "ベースライン時の年齢（歳）",
        "学歴（学校教育年数）",
        "ベースライン時の1日あたり喫煙本数",
        "ベースライン時の喫煙年数",
        "運動レベル（カテゴリ変数）",
        "身体活動の有無（カテゴリ変数）",
        "ベースライン時の体重（ポンド）",
        "1971年から1982年の体重変化（ポンド）"
    )
)
knitr::kable(var_nhefs_desc, caption = "nhefsデータセットの各変数の説明")

# アウトカムを1971年から1982年の体重変化（wt82_71）、主要な曝露変数を禁煙指標（qsmk）、その他の変数を共変量とし重回帰分析を行います。今回あえてtable1を作成していませんが、自主トレーニングとして、gtsummaryでぜひ作成してみてください。

# 重回帰分析は`lm()`で行います。 結果は以下で表示することが多いです。

# 1.  summary()で表示
# 2.  broom::tidy()でデータフレームとして表示
# 3.  tbl_regression()で表示

# broomパッケージは色々なモデルの結果をデータフレームとして表示でき、summary()では表示されない95%信頼区間やp値も表示できます。gtsummary:tbl_regression()はbroom::tidy()の結果を受け、いい感じの表を作ってくれます。

# 多重共線性の確認で使われるVIFは`add_vif()`で表示できます。

# add_vif()はGVIFとAdjusted GVIFを表示します。GVIFは多重共線性の指標で、10を超えると多重共線性の目安になります。Adjusted GVIFは共変量に連続変数とカテゴリ変数が混在する場合に使われるらしいですが、Adjusted GVIFを指定しているサイトや書籍をあまり見かけないので、ここではGVIFを表示します。

# `add_significance_stars()`はp値に星をつけて表示します。

# 重回帰分析 VIFも追加
model_nhefs <- lm(
    wt82_71 ~
        qsmk +
            sex +
            race +
            age +
            as.factor(education) +
            smokeintensity +
            smokeyrs +
            as.factor(exercise) +
            as.factor(active) +
            wt71,
    data = nhefs
)

summary_nhefs <- model_nhefs |>
    tbl_regression() |>
    add_significance_stars(hide_ci = FALSE, hide_se = TRUE, hide_p = FALSE) |>
    add_vif(statistic = "GVIF")

summary_nhefs

# もし因果推論として行う場合*table2の誤謬*を起こさないために曝露変数のみを表示することがあります。

# そういった場合は`include =` 引数を使い曝露変数のみを表示します。

# 重回帰分析の必要な行だけ抽出
summary_nhefs2 <-
    model_nhefs |>
    tbl_regression(include = qsmk)
summary_nhefs2

# tbl_regression()の後にplot()を使うと回帰係数のグラフを出力できます。

# グラフはggsave()を使うとグラフを保存できます。

# 回帰係数のグラフ
plot_nhefs <- summary_nhefs |>
    plot()
ggsave("plot_nhefs.png", plot_nhefs, width = 10, height = 8)

# その気になれば単変量解析と多変量解析の結果を並べることも可能です。

# まず各列に対して単変量解析を行い、その結果をデータフレームとして保存します。

# 単変量解析の結果をデータフレームとして保存
summary_nhefs_uv <-
    tbl_uvregression(
        nhefs,
        method = lm,
        y = wt82_71,
        include = qsmk
    )

tbl_merge(
    list(summary_nhefs_uv, summary_nhefs2),
    tab_spanner = c("**単変量解析**", "**多変量解析**")
)

## 2. ロジスティック回帰（pbcデータ）

# pbcデータを使いロジスティック回帰分析を行います。

# アウトカムのstatusを二値化してから解析します。

# statusが0を0, 1と2を1に変換します。

# ロジスティック回帰分析の結果はそのままだと回帰係数が表示されます。ただ解釈では回帰係数をオッズ比に変換することが多く、その設定が`exponentiate = TRUE`です。

# ロジスティック回帰分析
pbc$status_cat2 <- ifelse(pbc$status == 0, 0, 1)

model_logit <- glm(
    status_cat2 ~
        trt +
            sex +
            age +
            bili +
            albumin +
            copper +
            alk.phos +
            ast +
            trig +
            protime,
    data = pbc,
    family = binomial
)

summary_logit <- model_logit |>
    tbl_regression(exponentiate = TRUE) |>
    add_significance_stars(hide_ci = FALSE, hide_se = TRUE, hide_p = FALSE)

summary_logit

# グラフも同様にオッズ比に変換する方が解釈がしやすいです。

# ロジスティック回帰分析のグラフ
plot_logit <- summary_logit |>
    plot()
ggsave("plot_logit.png", plot_logit, width = 10, height = 8)

## 3. Cox比例ハザードモデル（pbcデータ）

# Coxモデルでも `exponentiate = TRUE` を指定することで、ハザード比（HR）を表示できます。これにより、解釈が直感的になり、報告にも適しています。

# Cox比例ハザードモデル
model_cox <- coxph(
    Surv(time, status) ~
        trt +
            sex +
            age +
            bili +
            albumin +
            copper +
            alk.phos +
            ast +
            trig +
            protime,
    data = pbc
)
summary_cox <- model_cox |>
    tbl_regression(exponentiate = TRUE) |>
    add_significance_stars(hide_ci = FALSE, hide_se = TRUE, hide_p = FALSE)

summary_cox

## 4. 傾向スコアマッチング（pbcデータ）

# 傾向スコアマッチングではマッチング前後でSDMを出力することがあります。

# pbcで傾向スコアマッチングを行います。

# 傾向スコアマッチング
theme_gtsummary_compact()
pbc_na_omit <- pbc |>
    na.omit()

pbc_match <- matchit(
    trt ~ sex + age + bili + albumin + copper + alk.phos + ast + trig + protime,
    data = pbc_na_omit,
    method = "nearest",
    link = "linear.logit",
    caliper = 0.2,
    std.caliper = TRUE,
    ratio = 1
)

match_data <- match.data(pbc_match)

summary_before_match <- pbc_na_omit |>
    tbl_summary(by = trt, include = include_vars) |>
    add_difference(everything() ~ "smd") |>
    modify_column_hide("conf.low")

summary_after_match <- match_data |>
    tbl_summary(by = trt, include = include_vars) |>
    add_difference(everything() ~ "smd") |>
    modify_column_hide("conf.low")

tbl_merge(
    list(summary_before_match, summary_after_match),
    tab_spanner = c("**マッチング前**", "**マッチング後**")
)

# マッチング後のロジスティック回帰分析を行います。ここではロバスト分散など考慮せずに行います。

model_match <- glm(status_cat2 ~ trt, data = match_data, family = binomial)
summary_match <- model_match |>
    tbl_regression(exponentiate = TRUE) |>
    add_significance_stars(hide_ci = FALSE, hide_se = TRUE, hide_p = FALSE)

summary_match
