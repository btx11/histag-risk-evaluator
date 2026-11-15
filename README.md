# de novo 酵素 Hisタグ付加リスク評価 & PyMOL 3方向バッチ作図ツール (`histag-risk-evaluator`)

ディレクトリ内の複数の **酵素 + リガンド（基質・補酵素の2種など）ドッキング立体構造 PDB ファイル** を一括解析し、**N末端とC末端のどちらに Hisタグ（$\text{His}_6$ 等）を付加するリスクが高いか** を構造生物学・タンパク質工学の先行研究（論文）に基づいて定量評価し、**PyMOL で色分けした 3方向（正面 0° + 90°回転×2種）の立体構造図** とともに **Excel (`.xlsx`)** に自動集約する Python スクリプトです。

Windows PC にインストールされた **デスクトップ版 PyMOL** と **WSL2 (Ubuntu)** の両方からシームレスに動作するよう設計されています。

---

## 主な機能

1. **PyMOL による自動色分け & 3方向バッチレンダリング（計3種の図 / PDB）**
   - **酵素（タンパク質鎖）**: 青色 (`blue`, Cartoon表示)
   - **リガンド（基質・補酵素の2種）**: 黄色 (`yellow`, Stick表示)
   - **N末端**: 赤色 (`red`, 末端残基 Sphere + 末端領域 Cartoon 表示)
   - **C末端**: 緑色 (`green`, 末端残基 Sphere + 末端領域 Cartoon 表示)
   - **3方向の画像を自動作成**:
     - `View 1`: `orient` による正面基準図 (0°)
     - `View 2`: Y軸周りに **90°回転** させた図 (`turn y, 90`)
     - `View 3`: X軸周りに **90°回転**（または Y軸周りにさらに90°回転、`--rot-mode` で選択可能）させた図
   - 各 PDB に対応する再現用 PyMOL スクリプト (`output_images/pml_scripts/*.pml`) も自動出力します。

2. **先行研究（論文）に基づく N末・C末の定量評価指標**
   - **リガンド（基質・補酵素）との距離 ($d_{\min}, d_{\text{C}\alpha\text{-centroid}}$)**
     - 全リガンドに対する最小重原子間距離 $d_{\min}$ ($\text{\AA}$) に加え、**基質 (Substrate)** と **補酵素 (Cofactor: NAD, FAD, HEM, SAM, ATP, PLP 等)** のそれぞれに対する個別距離も算出します。
     - $\text{His}_6$ タグの物理的コンター長（完全伸長時 $\sim 3.6\text{ \AA}/\text{残基} \times 6 \approx 21.6\text{ \AA}$、ランダムコイル実効到達半径 $\sim 12\text{–}15\text{ \AA}$）に基づき、活性部位ポケットへの立体障害・静電的干渉リスクを評価します（Carson et al., 2007; Freydank et al., 2008; Yeon et al., 2014）。
   - **溶媒露出表面積 (SASA) と相対溶媒露出度 (RSA)**
     - **Shrake & Rupley (1973)** の Rolling-ball 法（水分子プローブ半径 $1.4\text{ \AA}$、黄金螺旋 960 点）により各原子の絶対 $\text{SASA}$ ($\text{\AA}^2$) を計算します。
     - **Tien et al. (2013)** の理論的最大露出表面積 ($\text{MaxASA}$) により正規化した **相対溶媒露出度 ($\text{RSA} = \text{SASA} / \text{MaxASA}$)** を、最末端1残基 ($\text{RSA}_{\text{1res}}$) および末端3残基窓平均 ($\text{RSA}_{\text{win3}}$, Thornton & Sibanda, 1983) について算出します。
     - **Miller et al. (1987)** および **Rost & Sander (1994)** の基準に従い、埋没 ($\text{RSA} < 20\%$)、部分露出 ($20\% \le \text{RSA} < 50\%$)、露出 ($\text{RSA} \ge 50\%$) を判定します。
   - **末端主鎖ベクトルの配向角 ($\theta$)**
     - 末端から外側へ伸びる主鎖ベクトルと活性部位重心方向のなす角 $\theta$ ($0^\circ\text{–}180^\circ$) を算出し、付加したタグが活性部位側へ伸びやすいか溶媒側へ逃げるかを評価します。

3. **Excel (`.xlsx`) への自動集約（3シート構成）**
   - `Summary_&_3Views`: 各 PDB の推奨タグ位置、N末・C末の距離・SASA/RSA・リスクスコア・判定コメント、および **3方向の立体構造画像** を1行ずつ一覧表示。
   - `Detailed_Metrics`: 各鎖・各リガンド（基質・補酵素個別）との重原子最小距離・$\text{C}_\alpha$重心距離・B-factor/pLDDT・配向角・内訳スコアの詳細シート。
   - `Methodology_&_References`: 評価数式・閾値根拠および引用論文リスト（DOI付き）。

---

## 評価基準と引用文献 (Methodology & References)

### 1. 溶媒露出度 (SASA / RSA) の計算と埋没リスク判定
$$\text{RSA}_i = \min\left(1.0,\ \frac{\text{SASA}_i}{\text{MaxASA}_{\text{Tien2013}}(\text{AA}_i)}\right), \quad \text{RSA}_{\text{win3}} = \frac{1}{3}\sum_{k=1}^{3} \text{RSA}_k$$

| 相対溶媒露出度 ($\text{RSA}_{\text{win3}}$) | 分類 (Miller et al., 1987; Rost & Sander, 1994) | Hisタグ付加における構造・精製リスク |
| :--- | :--- | :--- |
| $\text{RSA}_{\text{win3}} < 0.20$ ($20\%$未満) | **Buried (埋没)** | **高リスク**: 末端が疎水性コアにパッキングされており、タグ直接付加によりフォールディング阻害・封入体化や、Ni-NTAレジンへの結合不全（タグの立体障害）を引き起こすリスクが高い。 |
| $0.20 \le \text{RSA}_{\text{win3}} < 0.50$ ($20\text{–}50\%$) | **Partially Exposed (部分露出)** | **中リスク**: 短いフレキシブルリンカー（例: $\text{GGGS}$ または $\text{(GGGS)}_2$）を挿入して溶媒側へ逃がす設計を推奨。 |
| $\text{RSA}_{\text{win3}} \ge 0.50$ ($50\%$以上) | **Exposed (溶媒露出)** | **低リスク**: 末端がタンパク質表面で十分に溶媒露出しており、Hisタグ付加に好適。 |

### 2. 基質・補酵素との距離 ($d_{\min}$) による活性部位干渉リスク判定
| リガンドまでの最小重原子距離 ($d_{\min}$) | 分類 (Carson et al., 2007; Freydank et al., 2008; Yeon et al., 2014) | 酵素活性・基質結合への干渉リスク |
| :--- | :--- | :--- |
| $d_{\min} < 15.0\text{ \AA}$ | **High Risk (高リスク)** | $\text{His}_6$ タグの一般的な揺らぎ到達圏内。基質侵入チャネルの立体閉塞、触媒ループの変形、イミダゾール環による補酵素・金属配位攪乱の危険が高い。 |
| $15.0\text{ \AA} \le d_{\min} < 22.0\text{ \AA}$ | **Moderate Risk (中リスク)** | $\text{His}_6$ タグの最大完全伸長コンター長（$\sim 21.6\text{ \AA}$）圏内。長すぎるリンカーは避け、配向角 $\theta$ を確認することを推奨。 |
| $d_{\min} \ge 22.0\text{ \AA}$ | **Low Risk (低リスク)** | $\text{His}_6$ タグの物理的到達圏外であり、活性部位への直接干渉リスクは極めて低い。 |

### 3. 総合リスクスコア ($\text{Total Risk Score} \in [0, 100]$)
$$\text{TotalRiskScore} = \max\Big(0.55 \cdot S_{\text{dist}}(d_{\min}) + 0.35 \cdot S_{\text{RSA}}(\text{RSA}_{\text{win3}}) + 0.10 \cdot S_{\text{orient}}(\theta),\ 0.85 \cdot \max(S_{\text{dist}}, S_{\text{RSA}})\Big)$$
（距離と溶媒露出度のいずれかが致命的である場合にボトルネック補正が働く設計となっています。）

### 引用論文一覧
1. **Shrake, A., & Rupley, J. A. (1973).** Environment and exposure to solvent of protein atoms. Lysozyme and insulin. *Journal of Molecular Biology*, 79(2), 351–371. DOI: [10.1016/0022-2836(73)90011-9](https://doi.org/10.1016/0022-2836(73)90011-9)
2. **Tien, M. Z., Meyer, A. G., Sydykova, D. K., Spielman, S. J., & Wilke, C. O. (2013).** Maximum allowed solvent accessibilities of residues in proteins. *PLOS ONE*, 8(11), e80635. DOI: [10.1371/journal.pone.0080635](https://doi.org/10.1371/journal.pone.0080635)
3. **Thornton, J. M., & Sibanda, B. L. (1983).** Amino and carboxy-terminal regions in globular proteins. *Journal of Molecular Biology*, 167(2), 443–460. DOI: [10.1016/S0022-2836(83)80344-1](https://doi.org/10.1016/S0022-2836(83)80344-1)
4. **Miller, S., Janin, J., Lesk, A. M., & Chothia, C. (1987).** Interior and surface of monomeric proteins. *Journal of Molecular Biology*, 196(3), 641–656. DOI: [10.1016/0022-2836(87)90038-6](https://doi.org/10.1016/0022-2836(87)90038-6)
5. **Rost, B., & Sander, C. (1994).** Conservation and prediction of solvent accessibility in protein families. *Proteins: Structure, Function, and Bioinformatics*, 20(3), 216–226. DOI: [10.1002/prot.340200303](https://doi.org/10.1002/prot.340200303)
6. **Carson, M., Johnson, D. H., McDonald, H., Brouillette, C., & DeLucas, L. J. (2007).** His-tag impact on structure. *Acta Crystallographica Section D: Biological Crystallography*, 63(3), 295–301. DOI: [10.1107/S0907444906052024](https://doi.org/10.1107/S0907444906052024)
7. **Freydank, A. C., Brandt, W., & Dräger, B. (2008).** Protein structure modeling indicates hexahistidine-tag interference with enzyme activity. *Proteins: Structure, Function, and Bioinformatics*, 72(1), 173–183. DOI: [10.1002/prot.21905](https://doi.org/10.1002/prot.21905)
8. **Yeon, Y. J., Park, H. J., & Yoo, Y. J. (2014).** Effect of His-tag location on the catalytic activity of 3-hydroxybutyrate dehydrogenase. *Biotechnology and Bioprocess Engineering*, 19(5), 798–802. DOI: [10.1007/s12257-014-0089-2](https://doi.org/10.1007/s12257-014-0089-2)

---

## セットアップ & 使い方 (Windows & WSL2 Ubuntu 両対応)

### 1. 依存パッケージのインストール（任意）
本スクリプトは **Python 標準ライブラリのみでも完全動作**（SASA/RSA計算および画像埋め込み `.xlsx` 出力の内蔵エンジンを搭載）しますが、`openpyxl` をインストールするとよりリッチなセル装飾で出力できます。

```bash
pip install -r requirements.txt
```

### 2. WSL2 (Ubuntu) から実行する場合（Windows デスクトップ版 PyMOL を自動呼び出し）
スクリプトは `/mnt/c/Users/<User>/AppData/Local/Schrodinger/PyMOL2/PyMOLWin.exe` 等にインストールされた **Windows デスクトップ版 PyMOL** を自動検出し、`wslpath` 変換および一時ステージングを介してヘッドレス (`-cq`) で 3方向の画像をバッチ生成します。

```bash
python3 evaluate_histag_risk.py \
  -i ./pdbs \
  -o ./histag_risk_evaluation.xlsx \
  --image-dir ./output_images
```

もし PyMOL のインストール先がカスタムパスの場合は `--pymol-cmd` で指定できます：
```bash
python3 evaluate_histag_risk.py \
  -i ./pdbs \
  -o ./histag_risk_evaluation.xlsx \
  --pymol-cmd "/mnt/c/Users/YourName/AppData/Local/Schrodinger/PyMOL2/PyMOLWin.exe"
```

### 3. Windows (PowerShell / コマンドプロンプト) から直接実行する場合
```powershell
py evaluate_histag_risk.py -i .\pdbs -o .\histag_risk_evaluation.xlsx --image-dir .\output_images
```

### 4. 基質 (Substrate) と補酵素 (Cofactor) の残基名を明示指定したい場合
デフォルトで標準的な補酵素（`NAD, NAP, FAD, FMN, HEM, SAM, ATP, PLP, COA` 等）とドッキング基質（`SUB, LIG, UNL` 等）を自動判別しますが、明示的に残基名を指定することも可能です：
```bash
python3 evaluate_histag_risk.py \
  -i ./pdbs \
  -o ./histag_risk_evaluation.xlsx \
  --substrate-resn SUB LIG \
  --cofactor-resn NAD FAD
```

### 5. サンプル PDB での動作確認
同梱の `generate_sample_pdbs.py` を実行すると、3種類のモデルケース（①N末が活性部位近傍かつ埋没、②C末が活性部位近傍、③両末端とも露出かつ遠위）のドッキングPDBを生成してすぐにテストできます：
```bash
python3 generate_sample_pdbs.py
python3 evaluate_histag_risk.py -i ./sample_pdbs -o ./sample_output/histag_risk_evaluation.xlsx --image-dir ./sample_output/images
```
