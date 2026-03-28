"""
pymna-pro2.py
Copyright (c) 2026 zyutama
Released under the MIT license
https://opensource.org/licenses/mit-license.php
"""
import numpy as np
import tkinter as tk
from tkinter import messagebox
import re
import matplotlib.pyplot as plt
from tkinter import simpledialog
import os

# --- 1. 単位変換エンジン ---
UNIT_DICT = {
    'g': 1e9, 'G': 1e9, 'meg': 1e6, 'Meg': 1e6, 'MEG': 1e6,
    'k': 1e3, 'K': 1e3, 'm': 1e-3, 'u': 1e-6, 'U': 1e-6,
    'n': 1e-9, 'N': 1e-9, 'p': 1e-12, 'P': 1e-12,
}

def parse_value(val_str):
    if not val_str or str(val_str).strip() in ['-', '', 'none', '0']:
        return 0.0
    val_str = str(val_str).strip()
    match = re.match(r"([+-]?\d*\.?\d+)([a-zA-Z]+)?", val_str)
    if not match: return 0.0
    num = float(match.group(1))
    unit = match.group(2)
    if unit and unit in UNIT_DICT:
        num *= UNIT_DICT[unit]
    return num

# --- 2. 修正節点解析 (MNA) エンジン ---
class MNASolver:
    def __init__(self, netlist):
        self.netlist = netlist
        self.nodes = set()
        for row in netlist:
            if len(row) < 4: continue
            self.nodes.add(str(row[2]))
            self.nodes.add(str(row[3]))
            # For Voltage Controlled Voltage Source (E) and Current Controlled Current Source (F)
            # need to add controlling nodes to node map, though 'F' refers to a V-source name.
            if row[0].upper() == 'E' and len(row) >= 6:
                self.nodes.add(str(row[4])) # Controlling positive node for E source
                self.nodes.add(str(row[5])) # Controlling negative node for E source
            # F type (CCCS) format: F Name N+ N- VCTRL GAIN
            # VCTRL is a voltage source name, not a node itself, so no need to add row[4] as node

        self.nodes.discard('0')
        self.node_map = {name: i for i, name in enumerate(sorted(list(self.nodes)))}
        self.num_n = len(self.node_map)
        # Include 'V' and 'E' type sources for MNA extra variables (their currents)
        # 'F' type sources need to reference currents of existing 'V' or 'E' sources.
        self.v_sources = [r for r in netlist if r[0].upper() in ['V', 'E']] # 'F' type removed from here
        self.num_v = len(self.v_sources)
        self.v_map = {r[1]: i for i, r in enumerate(self.v_sources)}
        self.dim = self.num_n + self.num_v

    def solve(self, params):
        A = np.zeros((self.dim, self.dim))
        Z = np.zeros(self.dim)
        for row in self.netlist:
            rtype, name = row[0].upper(), row[1]
            # The `val` from params will be used for component value or gain.
            # For resistors, current sources, voltage sources, and E-sources, it's directly used.
            # For F-sources, `val` will be the gain.
            val = params.get(name, 0.0)
            n_pos, n_neg = self.node_map.get(str(row[2]), -1), self.node_map.get(str(row[3]), -1)

            if rtype == 'R':
                g = 1.0 / val if val != 0 else 1e12
                if n_pos != -1: A[n_pos, n_pos] += g
                if n_neg != -1: A[n_neg, n_neg] += g
                if n_pos != -1 and n_neg != -1: A[n_pos, n_neg] -= g; A[n_neg, n_pos] -= g
            elif rtype == 'I':
                # Current source flowing from n_neg to n_pos
                if n_pos != -1: Z[n_pos] -= val
                if n_neg != -1: Z[n_neg] += val
            elif rtype == 'V':
                # Voltage source from n_pos to n_neg
                v_idx = self.num_n + self.v_map[name]
                if n_pos != -1: A[n_pos, v_idx] += 1; A[v_idx, n_pos] += 1
                if n_neg != -1: A[n_neg, v_idx] -= 1; A[v_idx, n_neg] -= 1
                Z[v_idx] = val
            elif rtype == 'E':
                # Voltage Controlled Voltage Source (VCVS)
                # E Name N+ N- NI+ NI- GAIN
                v_idx = self.num_n + self.v_map[name]
                ni_pos, ni_neg = self.node_map.get(str(row[4]), -1), self.node_map.get(str(row[5]), -1)
                
                if n_pos != -1: A[n_pos, v_idx] += 1; A[v_idx, n_pos] += 1
                if n_neg != -1: A[n_neg, v_idx] -= 1; A[v_idx, n_neg] -= 1
                
                # Equation for VCVS: V(N+,N-) - GAIN * V(NI+,NI-) = 0
                # -> V_v_idx - GAIN * (V_ni_pos - V_ni_neg) = 0
                # In MNA, this is implemented by modifying the row corresponding to V_v_idx
                if ni_pos != -1: A[v_idx, ni_pos] -= val # -= GAIN * V(NI+)
                if ni_neg != -1: A[v_idx, ni_neg] += val # += GAIN * V(NI-)
            elif rtype == 'F':
                # Current Controlled Current Source (CCCS)
                # F Name N+ N- VCTRL GAIN
                # Current of GAIN * I(VCTRL) flows from N+ to N-
                vctrl_name = row[4] # Name of the controlling voltage source
                gain = val # GAIN is the component's value for F type

                if vctrl_name not in self.v_map: # Check if controlling V source is defined
                    raise ValueError(f"Controlling voltage source '{vctrl_name}' for CCCS '{name}' not found.")

                # Get the index corresponding to the current of the controlling voltage source I(VCTRL)
                vctrl_idx = self.num_n + self.v_map[vctrl_name]

                # Add gain to the MNA matrix for the current flowing from N+ to N-
                if n_pos != -1: A[n_pos, vctrl_idx] += gain
                if n_neg != -1: A[n_neg, vctrl_idx] -= gain

        try:
            # Debugging: Print A and Z before solving
            print("--- MNA Matrix A ---")
            print(A)
            print("--- MNA Vector Z ---")
            print(Z)
            print(f"Node Map: {self.node_map}")
            print(f"Voltage Source Map: {self.v_map}")

            res = np.linalg.solve(A, Z)
            volts = {n: res[self.node_map[n]] for n in self.node_map}
            volts['0'] = 0.0
            return volts
        except Exception as e:
            print(f"Error solving MNA: {e}") # For debugging - RE-ENABLED THIS LINE
            return {n: 0.0 for n in self.node_map} | {'0': 0.0}


# --- 3. 言語リソース ---
TEXTS = {
    "jp": {
        "guide_title": "📌 入力ガイド",
        "guide_text": "【基本書式】\nR/V/I  Name  n+  n-  Typ  [Tol/Range]\n\n● 許容誤差(%)で指定する場合\n  R  R1  1  0  1k  1  (±1%)\n\n● 最小/最大を実値で指定する場合\n  R  R1  1  0  1k  0.95k/1.1k\n  (非対称なマージン設定が可能)\n\n【オペアンプ(E)書式】\nE  Name  o+  o-  i+  i-  Gain\n例: E OP1 3 0 1 2 100k\n\n【電流制御電流源(F)書式】\nF  Name  o+  o-  VCTRL  Gain  [Tol/Range]\n例: F F1 3 0 Vx 2 1.8/2.2\n（VCTRLは制御する電流源名）\n\n【単位】\nk, meg, m, u, n, p 対応\n\n※Excelからコピー＆ペースト可",
        "preset_lbl": "プリセット回路",
        "load_btn": "ロード",
        "lt_check": "LTspice保存",
        "iter_check": "反復収束法",
        "run_btn": "解析実行",
        "res_header": "--- Node Voltages ---\nNode\tTyp[V]\tMin[V]\tMax[V]\n",
        "presets": {
            "div": "R-2Rラダー回路",
            "super": "重ね合わせ回路",
            "inv": "反転増幅回路",
            "noninv": "非反転増幅回路"
        },
        "toggle_lang": "English Interface"
    },
    "en": {
        "guide_title": "📌 Input Guide",
        "guide_text": "[Basic Format]\nR/V/I  Name  n+  n-  Typ  [Tol/Range]\n\n* Tolerance(%):\n  R  R1  1  0  1k  1  (±1%)\n\n* Min/Max Values:\n  R  R1  1  0  1k  0.95k/1.1k\n  (Asymmetric margins allowed)\n\n[Op-Amp (E) Format]\nE  Name  o+  o-  i+  i-  Gain\nEx: E OP1 3 0 1 2 100k\n\n[CCCS (F) Format]\nF  Name  o+  o-  VCTRL  Gain  [Tol/Range]\nEx: F F1 3 0 Vx 2 1.8/2.2\n(VCTRL is controlling source name)\n\n[Units]\nk, meg, m, u, n, p supported\n\n* Copy & Paste from Excel available",
        "preset_lbl": "Preset Circuits",
        "load_btn": "Load",
        "lt_check": "Save LTspice",
        "iter_check": "Iterative Refinement",
        "run_btn": "RUN ANALYSIS",
        "res_header": "--- Node Voltages ---\nNode\tTyp[V]\tMin[V]\tMax[V]\n",
        "presets": {
            "div": "R-2R Ladder",
            "super": "Superposition",
            "inv": "Inverting Amp",
            "noninv": "Non-Inverting Amp"
        },
        "toggle_lang": "日本語インターフェース"
    }
}

PRESET_DATA_RAW = {
    "div": {"netlist": "V\tVin\t1\t0\t10\nR\tR1\t1\t2\t1k\t1\nR\tR2\t1\t0\t2k\t1\nR\tR3\t2\t3\t1k\t1\nR\tR4\t2\t0\t2k\t1\nR\tR5\t3\t4\t1k\t1\nR\tR6\t3\t0\t2k\t1\nR\tR7\t4\t5\t1k\t1\nR\tR8\t4\t0\t2k\t1\nR\tR9\t5\t0\t1k\t1", "ascii": "Vin(1)o--+--R1(1k)--+--R3(1k)--+--R5(1k)--+--R7(1k)--o(5)\n         |          |          |          |          |\n       R2(2k)     R4(2k)     R6(2k)     R8(2k)     R9(1k)\n         |          |          |          |          |\n        GND        GND        GND        GND        GND"},
    "super": {"netlist": "V\tV1\t1\t0\t10\nV\tV2\t2\t0\t5\nR\tR1\t1\t3\t1k\t2\nR\tR2\t2\t3\t2k\t2\nR\tR3\t3\t0\t1k\t2", "ascii": "V1(10V) o---R1---+\n                 V2(5V) o---R2---+\n                              |\n                              R3\n                              |\n                              0"},
    "inv": {"netlist": "V\tVin\t1\t0\t1\nR\tRin\t1\t2\t10k\t1\nR\tRf\t2\t3\t100k\t1\nE\tOP1\t3\t0\t2\t0\t100k", "ascii": "Vin--Rin--(-)\n            |       |\n            +--Rf---(+)/0"},
    "noninv": {"netlist": "V\tVin\t1\t0\t1\nR\tR1\t3\t2\t10k\t1\nR\tR2\t2\t0\t100k\t1\nE\tOP1\t3\t0\t1\t2\t100k", "ascii": "Vin--+\n      |\n      R1\n      |\n      +--R2--0\n      |\n     (+)"}
}

# --- 3. メインGUI ---
class PyMNAProApp(tk.Tk):
    def save_ltspice_netlist(self, netlist, comps):
        """LTspiceでそのまま動くテキストネットリスト（.cir/.net）を出力"""
        try:
            lines = []
            # 部品記号変換辞書
            symbol_map = {'R': 'R', 'V': 'V', 'I': 'I', 'E': 'E', 'F': 'F'}
            for row in netlist:
                rtype = row[0].upper()
                name = row[1]
                n1 = row[2]
                n2 = row[3]
                if rtype in ['R', 'V', 'I']:
                    value = comps[name]['raw']
                    lines.append(f"{symbol_map[rtype]}{name} {n1} {n2} {value}")
                elif rtype == 'E':
                    # E Name N+ N- NI+ NI- GAIN
                    ni_pos = row[4]
                    ni_neg = row[5]
                    gain = comps[name]['raw']
                    lines.append(f"E{name} {n1} {n2} {ni_pos} {ni_neg} {gain}")
                elif rtype == 'F':
                    # F Name N+ N- VCTRL GAIN
                    vctrl = row[4]
                    gain = comps[name]['raw']
                    lines.append(f"F{name} {n1} {n2} V{vctrl} {gain}")
                # 他の素子は必要に応じて追加
            # 解析制御行（例: 適当な.tran）
            lines.append('.tran 0 1 0 0.01')
            lines.append('.end')
            with open("output_circuit.cir", "w", encoding="utf-8") as f:
                f.write("* LTspice Netlist Export\n")
                for l in lines:
                    f.write(l + "\n")
            self.output_text.insert("end", "\n[System] LTspice netlist (.cir) saved.")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def __init__(self):
        super().__init__()
        self.lang = "en" # Default to English
        self.title("pymna-pro")
        self.geometry("1400x950")
        self.configure(bg="#23272e")

        self.grid_columnconfigure(0, weight=0)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        # サイドバー
        self.sidebar = tk.Frame(self, width=320, bg="#23272e", bd=0, highlightthickness=0)
        self.sidebar.grid(row=0, column=0, sticky="nsew")
        
        # 言語切り替えボタン
        self.btn_lang = tk.Button(self.sidebar, text=TEXTS[self.lang]["toggle_lang"], command=self.toggle_language, font=("Segoe UI", 9), bg="#444", fg="white", relief="flat")
        self.btn_lang.pack(anchor="ne", padx=5, pady=5)

        self.lbl_guide_title = tk.Label(self.sidebar, text=TEXTS[self.lang]["guide_title"], font=("Segoe UI", 20, "bold"), fg="#00bfff", bg="#23272e", anchor="w")
        self.lbl_guide_title.pack(pady=(10,8), padx=16, anchor="w")

        self.lbl_guide = tk.Label(self.sidebar, text=TEXTS[self.lang]["guide_text"], justify="left", font=("Consolas", 12), anchor="w", fg="#e0e0e0", bg="#23272e")
        self.lbl_guide.pack(padx=16, pady=(0,8), fill="x")

        # プリセットUIコンテナ
        self.preset_frame = tk.Frame(self.sidebar, bg="#23272e")
        self.preset_frame.pack(fill="x", padx=16, pady=(0,8))
        self.preset_var = tk.StringVar(value="div") # IDで管理
        self.update_preset_ui() # 初回構築

        self.ascii_label = tk.Label(self.sidebar, text="", justify="left", anchor="nw", font=("Consolas", 9), fg="#b0ffb0", bg="#23272e")
        self.ascii_label.pack(fill="x", padx=16, pady=(0,14))

        # メインフレーム
        self.main_frame = tk.Frame(self, bg="#23272e")
        self.main_frame.grid(row=0, column=1, sticky="nsew", padx=0, pady=0)
        self.main_frame.grid_columnconfigure((0, 1), weight=1)
        self.main_frame.grid_rowconfigure(1, weight=1)

        self.input_text = tk.Text(self.main_frame, font=("Consolas", 13), borderwidth=0, relief="flat", bg="#181c22", fg="#e0e0e0", insertbackground="#00bfff", selectbackground="#2d3748")
        self.input_text.grid(row=1, column=0, sticky="nsew", padx=(12,6), pady=12)
        self.load_preset()

        self.output_text = tk.Text(self.main_frame, font=("Consolas", 13), borderwidth=0, relief="flat", bg="#101216", fg="#00ff99", insertbackground="#00bfff", selectbackground="#2d3748")
        self.output_text.grid(row=1, column=1, sticky="nsew", padx=(6,12), pady=12)

        self.ctrl = tk.Frame(self.main_frame, bg="#23272e")
        self.ctrl.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(0,10))
        self.iter_var = tk.BooleanVar(value=True)
        self.lt_var = tk.BooleanVar(value=True)
        self.chk_iter = tk.Checkbutton(self.ctrl, text=TEXTS[self.lang]["iter_check"], variable=self.iter_var, font=("Segoe UI", 11), fg="#e0e0e0", bg="#23272e", activebackground="#23272e", selectcolor="#23272e")
        self.chk_iter.pack(side="left", padx=20)
        self.chk_lt = tk.Checkbutton(self.ctrl, text=TEXTS[self.lang]["lt_check"], variable=self.lt_var, font=("Segoe UI", 11), fg="#e0e0e0", bg="#23272e", activebackground="#23272e", selectcolor="#23272e")
        self.chk_lt.pack(side="left", padx=20)
        self.btn_run = tk.Button(self.ctrl, text=TEXTS[self.lang]["run_btn"], command=self.execute, font=("Segoe UI", 12, "bold"), bg="#00bfff", fg="#23272e", activebackground="#0099cc", activeforeground="#ffffff", relief="flat", bd=0, width=18, height=2, cursor="hand2")
        self.btn_run.pack(side="right", padx=20)

        self.btn_tornado = tk.Button(self.ctrl, text="Tornado Chart", command=self.generate_tornado_chart, font=("Segoe UI", 11), bg="#555", fg="white", relief="flat", state="disabled")
        self.btn_tornado.pack(side="right", padx=(5, 0))

        self.solver = None
        self.comps = None
        self.v_typ = None

    def toggle_language(self):
        self.lang = "jp" if self.lang == "en" else "en"
        t = TEXTS[self.lang]
        self.btn_lang.config(text=t["toggle_lang"])
        self.lbl_guide_title.config(text=t["guide_title"])
        self.lbl_guide.config(text=t["guide_text"])
        self.chk_iter.config(text=t["iter_check"])
        self.chk_lt.config(text=t["lt_check"])
        self.btn_run.config(text=t["run_btn"])
        self.update_preset_ui()

    def update_preset_ui(self):
        for widget in self.preset_frame.winfo_children():
            widget.destroy()
        
        t = TEXTS[self.lang]
        tk.Label(self.preset_frame, text=t["preset_lbl"], font=("Segoe UI", 11, "bold"), fg="#ffffff", bg="#23272e").pack(anchor="w")
        
        display_map = t["presets"]
        display_names = list(display_map.values())
        current_id = self.preset_var.get()
        current_name = display_map.get(current_id, display_names[0])
        
        self.disp_var = tk.StringVar(value=current_name)
        
        def on_change(selected_name):
            for kid, vname in display_map.items():
                if vname == selected_name:
                    self.preset_var.set(kid)
                    break
        
        om = tk.OptionMenu(self.preset_frame, self.disp_var, *display_names, command=on_change)
        om.config(bg="#444", fg="white", highlightthickness=0, relief="flat")
        om["menu"].config(bg="#444", fg="white")
        om.pack(fill="x", pady=2)
        
        tk.Button(self.preset_frame, text=t["load_btn"], command=self.load_preset, font=("Segoe UI", 10), bg="#00bfff", fg="#23272e").pack(fill="x", pady=4)

    def load_preset(self):
        key = self.preset_var.get()
        item = PRESET_DATA_RAW.get(key)
        if not item:
            return
        self.input_text.delete("1.0", "end")
        self.input_text.insert("1.0", item["netlist"])
        self.ascii_label.config(text=item["ascii"])

    def generate_tornado_chart(self):
        if not self.solver or not self.v_typ:
            messagebox.showinfo("Info", "Please run a successful analysis first.")
            return

        target_node = simpledialog.askstring("Input", "Enter the node to analyze for the Tornado Chart:", parent=self)
        if not target_node or target_node not in self.v_typ:
            messagebox.showerror("Error", f"Invalid node '{target_node}'. Please enter a valid node from the analysis results.")
            return

        sensitivities = []
        baseline_v = self.v_typ[target_node]
        typ_params = {n: d['typ'] for n, d in self.comps.items()}

        # Calculate sensitivity for each component
        for name, comp_data in self.comps.items():
            if comp_data['min'] == comp_data['max']:
                continue

            # Calculate voltage at min value
            params_min = typ_params.copy()
            params_min[name] = comp_data['min']
            v_low = self.solver.solve(params_min)[target_node]

            # Calculate voltage at max value
            params_max = typ_params.copy()
            params_max[name] = comp_data['max']
            v_high = self.solver.solve(params_max)[target_node]
            
            delta_low = v_low - baseline_v
            delta_high = v_high - baseline_v
            swing = abs(delta_high - delta_low)
            
            sensitivities.append({'name': name, 'low': delta_low, 'high': delta_high, 'swing': swing})

        if not sensitivities:
            messagebox.showinfo("Info", "No components with tolerances found to generate a chart.")
            return

        sensitivities.sort(key=lambda x: x['swing'], reverse=True)

        comp_names = [s['name'] for s in sensitivities]
        low_values = np.array([min(s['low'], s['high']) for s in sensitivities])
        high_values = np.array([max(s['low'], s['high']) for s in sensitivities])

        fig, ax = plt.subplots(figsize=(10, max(6, len(comp_names) * 0.5)))
        bar_widths = high_values - low_values
        ax.barh(comp_names, bar_widths, left=low_values, color='skyblue', edgecolor='black', zorder=3)
        ax.axvline(0, color='red', linestyle='--', linewidth=1, zorder=4)

        ax.set_xlabel(f'Voltage Deviation from Typical ({baseline_v:.4f} V)')
        ax.set_ylabel('Component')
        ax.set_title(f'Tornado Chart: Sensitivity of Node "{target_node}" Voltage', fontsize=14, fontweight='bold')
        ax.grid(axis='x', linestyle=':', alpha=0.7, zorder=0)
        ax.invert_yaxis()
        fig.tight_layout()
        plt.show()

    def execute(self):
        try:
            raw = self.input_text.get("1.0", "end-1c").strip()
            if not raw: return

            # Improved netlist parsing to ignore comments
            lines = []
            for l in raw.split('\n'):
                l_strip = l.strip()
                if not l_strip or l_strip.startswith(';') or l_strip.startswith('*'): # Ignore empty or comment-only lines
                    continue
                # Remove inline comments starting with ';'
                if ';' in l_strip:
                    l_strip = l_strip.split(';', 1)[0].strip()
                if l_strip: # Only add if not empty after comment removal
                    lines.append(l_strip)

            netlist = [re.split(r'\t|\s+', l) for l in lines]
            solver = MNASolver(netlist)

            comps = {}
            for row in netlist:
                rtype, name = row[0].upper(), row[1]
                val = 0.0
                raw_val_str = ''
                v_min, v_max = 0.0, 0.0

                if rtype in ['R', 'V', 'I']:
                    val = parse_value(row[4])
                    raw_val_str = row[4]
                    if len(row) > 5:
                        tol_str = str(row[5]).strip()
                        if '/' in tol_str:
                            parts = tol_str.split('/')
                            v_min, v_max = parse_value(parts[0]), parse_value(parts[1])
                        else:
                            tol = parse_value(tol_str)
                            v_min, v_max = val*(1-tol/100), val*(1+tol/100)
                    else:
                        v_min, v_max = val, val
                elif rtype == 'E': # E OP1 3 0 1 2 100k (Gain is row[6])
                    val = parse_value(row[6]) if len(row) > 6 else parse_value(row[4])
                    raw_val_str = row[6] if len(row) > 6 else row[4]
                    v_min, v_max = val, val # Assuming fixed gain for E type for simplicity
                elif rtype == 'F': # F F1 N+ N- VCTRL GAIN (Gain is row[5])
                    val = parse_value(row[5])
                    raw_val_str = row[5]
                    v_min, v_max = val, val # Assuming fixed gain for F type for simplicity

                comps[name] = {'typ': val, 'min': v_min, 'max': v_max, 'type': rtype, 'raw': raw_val_str}

            v_typ = solver.solve({n: d['typ'] for n, d in comps.items()})
            report = TEXTS[self.lang]["res_header"] + "-"*45 + "\n"
            log = "\n--- Parameter Assignments ---\n"

            # 自然順（1, 2, ..., 11）でソート
            natural_sort_key = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]

            for target in sorted(solver.node_map.keys(), key=natural_sort_key):
                # 1. 独立した初期状態の決定 (感度解析)
                max_states = {}
                min_states = {}
                for n, d in comps.items():
                    if d['min'] == d['max']:
                        max_states[n] = min_states[n] = 'typ'
                    else:
                        tp = {k: v['typ'] for k, v in comps.items()}
                        tp[n] = d['max']
                        v_shifted = solver.solve(tp)[target]
                        if v_shifted > v_typ[target]:
                            max_states[n] = 'max'; min_states[n] = 'min'
                        else:
                            max_states[n] = 'min'; min_states[n] = 'max'

                if self.iter_var.get():
                    # --- MAX独立探索 ---
                    for _ in range(15):
                        changed = False
                        for n in comps:
                            if comps[n]['min'] == comps[n]['max']: continue
                            p_curr = {k: (comps[k]['max'] if max_states[k]=='max' else (comps[k]['min'] if max_states[k]=='min' else comps[k]['typ'])) for k in comps}
                            v_b = solver.solve(p_curr)[target]
                            old = max_states[n]; max_states[n] = 'min' if old=='max' else 'max' # Flip for trial
                            p_new = {k: (comps[k]['max'] if max_states[k]=='max' else (comps[k]['min'] if max_states[k]=='min' else comps[k]['typ'])) for k in comps}
                            if solver.solve(p_new)[target] > v_b: changed = True
                            else: max_states[n] = old # Revert if not better
                        if not changed: break

                    # --- MIN独立探索 ---
                    for _ in range(15):
                        changed = False
                        for n in comps:
                            if comps[n]['min'] == comps[n]['max']: continue
                            p_curr = {k: (comps[k]['max'] if min_states[k]=='max' else (comps[k]['min'] if min_states[k]=='min' else comps[k]['typ'])) for k in comps}
                            v_b = solver.solve(p_curr)[target]
                            old = min_states[n]; min_states[n] = 'min' if old=='max' else 'max'
                            p_new = {k: (comps[k]['max'] if min_states[k]=='max' else (comps[k]['min'] if min_states[k]=='min' else comps[k]['typ'])) for k in comps}
                            if solver.solve(p_new)[target] < v_b: changed = True # 減少を目指す
                            else: min_states[n] = old
                        if not changed: break

                p_max = {k: (comps[k]['max'] if max_states[k]=='max' else (comps[k]['min'] if max_states[k]=='min' else comps[k]['typ'])) for k in comps}
                v_max = solver.solve(p_max)[target]
                p_min = {k: (comps[k]['max'] if min_states[k]=='max' else (comps[k]['min'] if min_states[k]=='min' else comps[k]['typ'])) for k in comps}
                v_min = solver.solve(p_min)[target]

                report += f"{target}\t{v_typ[target]:.4f}\t{v_min:.4f}\t{v_max:.4f}\n"
                log += f"Node {target} (MAX selection): " + ", ".join([f"{k}:{max_states[k]}" for k in max_states if max_states[k]!='typ']) + "\n"
                log += f"Node {target} (MIN selection): " + ", ".join([f"{k}:{min_states[k]}" for k in min_states if min_states[k]!='typ']) + "\n"

            self.output_text.delete("1.0", "end")
            self.output_text.insert("1.0", report + log)

            self.solver = solver
            self.comps = comps
            self.v_typ = v_typ
            self.btn_tornado.config(state="normal", bg="#4CAF50")

            if self.lt_var.get():
                self.save_asc(netlist, comps)
                self.save_ltspice_netlist(netlist, comps)

        except Exception as e:
            self.btn_tornado.config(state="disabled", bg="#555")
            messagebox.showerror("Error", str(e))

    def save_asc(self, netlist, comps):
        try:
            with open("output_circuit.asc", "w", encoding="utf-8") as f:
                f.write("Version 4\nSHEET 1 880 680\n")
                x, y = 112, 96
                # LTspice標準のシンボル名に変換する辞書
                sym_map = {'r': 'res', 'v': 'voltage', 'i': 'current', 'c': 'cap', 'l': 'ind'}
                for row in netlist:
                    rtype, name = row[0].lower(), row[1]
                    sym = sym_map.get(rtype, rtype) # 辞書になければ元の文字(e, fなど)を使用
                    f.write(f"SYMBOL {sym} {x} {y} R0\nSYMATTR InstName {name}\nSYMATTR Value {comps[name]['raw']}\n")
                    y += 112
            self.output_text.insert("end", "\n[System] LTspice file saved.")
        except Exception as e:
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    # Removed the hardcoded netlist_str and reverted to GUI input for customtkinter version
    PyMNAProApp().mainloop()
