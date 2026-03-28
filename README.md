# pymna-pro-v2

**Automated Circuit WCA Solver: Finding exact worst-case nodes using MNA and iterative refinement.**
**自動回路WCAソルバー：MNAと反復改良を用いて正確な最悪ケースノードを探索**

## 概要 / Overview
本プロジェクトは、電子回路の信頼性設計において重要な**最悪値解析（WCA: Worst Case Analysis）**を効率化するためのPythonベースのシミュレータです。
従来の表計算ソフトのソルバー機能で課題となっていた収束精度や精度の問題を解決し、理論的な最悪値を確実に算出することを目的に開発されました。

This project is a Python-based circuit simulator designed to streamline **Worst Case Analysis (WCA)**, which is critical in reliable electronics design. It was developed to solve convergence and precision issues often encountered with general-purpose spreadsheet solvers, ensuring the calculation of theoretical worst-case values.

### 特徴 / Features
- **高精度な解法 / High Precision**: 修正節点解析（MNA）エンジンを搭載し、回路方程式を直接解くことで正確な節点電位を算出します。
  - Uses a Modified Nodal Analysis (MNA) engine to calculate exact node potentials by solving circuit equations directly.
- **自動最悪値探索 / Automated WCA**: 各部品の公差に基づき、特定のノード電圧を最大・最小にする組み合わせを反復収束アルゴリズムで自動的に見つけ出します。
  - Automatically identifies component value combinations that maximize or minimize specific node voltages based on tolerances using an iterative refinement algorithm.
- **感度分析の可視化 / Sensitivity Visualization**: トルネードチャート生成機能を備え、どの部品が回路特性に最も影響を与えているかを一目で把握できます。
  - Features Tornado Chart generation to visualize which components have the most significant impact on circuit characteristics.
- **実用的な連携 / Practical Integration**: ネットリストは汎用的な形式を採用しており、LTspice等へのエクスポート（.asc / .cir）も可能です。
  - Supports standard netlist formats and allows exporting to LTspice (.asc / .cir).

## 主な機能 / Key Functions
1. **MNAソルバー**: 抵抗、電圧源、電流源、連動電源（E, Fソース）に対応。
2. **WCA解析**: 各ノードの Typ / Min / Max 電圧を算出。
3. **トルネードチャート**: 特定ノードに対する各部品の寄与度をグラフ化。
4. **GUIインターフェース**: 直感的な操作と、表計算ソフトからのコピー＆ペーストによる入力に対応。

## セットアップ / Setup
```bash
pip install numpy matplotlib
