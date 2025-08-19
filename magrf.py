import tkinter as tk
import ttkbootstrap as tb
from ttkbootstrap.dialogs import Messagebox
import numpy as np
import cmath
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Circle
import matplotlib.patches as patches

# =============================================================================
# BLOCO DAS FUNÇÕES DE CÁLCULO (LÓGICA DO BACK-END)
# =============================================================================

C0 = 299792458  # Velocidade da luz no vácuo (m/s)

def dimensionar_microstrip(z0_ohms, h_metros, er):
    """
    Calcula a largura (W) de uma linha de microfita para um dado Z0, h e εr.
    """
    if z0_ohms <= 0 or h_metros <= 0 or er <= 1:
        return None
    A = (z0_ohms / 60) * np.sqrt((er + 1) / 2) + ((er - 1) / (er + 1)) * (0.23 + 0.11 / er)
    B = (377 * np.pi) / (2 * z0_ohms * np.sqrt(er))
    wh_larga_trial = (2 / np.pi) * (B - 1 - np.log(2 * B - 1) + ((er - 1) / (2 * er)) * (np.log(B - 1) + 0.39 - (0.61 / er)))
    if wh_larga_trial > 2:
        wh = wh_larga_trial
    else:
        wh = (8 * np.exp(A)) / (np.exp(2 * A) - 2)
    return wh * h_metros

def calcular_stub_paralelo_curto(frequencia_hz, z0_ohms, zl_ohms):
    """
    Calcula as dimensões (d, l) para um stub paralelo em curto-circuito.
    """
    if np.isclose(zl_ohms, z0_ohms):
        return {'solucoes': None, 'mensagem': "A impedância da carga já está casada (ZL = Z0)."}

    lambda_0 = C0 / frequencia_hz
    beta = (2 * np.pi) / lambda_0
    zl_norm = zl_ohms / z0_ohms
    yl_norm = 1 / zl_norm
    g = yl_norm.real
    b = yl_norm.imag
    
    A_quad = g**2 + b**2 - g
    B_quad = -2 * b
    C_quad = 1 - g

    if np.isclose(g, 1.0):
        t_solucoes = [0, 2 / b] if not np.isclose(b, 0) else [0]
    else:
        if np.isclose(A_quad, 0):
            if np.isclose(B_quad, 0):
                return {'solucoes': None, 'mensagem': "Não foi possível encontrar uma solução matemática."}
            t_solucoes = [-C_quad / B_quad]
        else:
            t_solucoes = np.roots([A_quad, B_quad, C_quad])

    solucoes_possiveis = []
    for t in t_solucoes:
        if np.iscomplex(t): continue
        d = np.arctan(t) / beta
        if d < 0: d += lambda_0 / 2
        
        yd_norm = (yl_norm + 1j * t) / (1 + 1j * yl_norm * t)
        b_prime = yd_norm.imag
        
        if np.isclose(b_prime, 0):
            l = lambda_0 / 4
        else:
            l = np.arctan(1 / b_prime) / beta
        
        if l < 0: l += lambda_0 / 2
        solucoes_possiveis.append({'d': d, 'l': l})

    if not solucoes_possiveis:
        return {'solucoes': None, 'mensagem': "Nenhuma solução real foi encontrada."}
    
    solucoes_possiveis.sort(key=lambda sol: sol['d'])
    
    return {
        'frequencia_hz': frequencia_hz, 'z0_ohms': z0_ohms, 'zl_ohms': zl_ohms,
        'lambda': lambda_0, 'yl_norm': yl_norm, 
        'solucoes': solucoes_possiveis[0]
    }

def calcular_s11_vs_frequencia(freq_min_ghz, freq_max_ghz, freq_step_ghz, z0_ohms, zl_ohms, solucao_stub):
    """
    Calcula o parâmetro S11 (perda de retorno) em função da frequência.
    """
    if not solucao_stub or 'd' not in solucao_stub or 'l' not in solucao_stub:
        return None
        
    d = solucao_stub['d']
    l = solucao_stub['l']
    freq_projeto_hz = solucao_stub.get('frequencia_hz', 0)
    freq_range_ghz = np.arange(freq_min_ghz, freq_max_ghz + freq_step_ghz, freq_step_ghz)
    s11_db = np.zeros(len(freq_range_ghz))
    
    for i, freq_ghz in enumerate(freq_range_ghz):
        freq_hz = freq_ghz * 1e9
        lambda_atual = C0 / freq_hz
        beta = 2 * np.pi / lambda_atual
        y_stub = (-1j / np.tan(beta * l)) / z0_ohms if np.tan(beta * l) != 0 else -1j * float('inf')
        z_in_d = z0_ohms * (zl_ohms + 1j * z0_ohms * np.tan(beta * d)) / (z0_ohms + 1j * zl_ohms * np.tan(beta * d))
        y_in_d = 1 / z_in_d
        y_total = y_in_d + y_stub
        z_total = 1 / y_total
        gamma_in = (z_total - z0_ohms) / (z_total + z0_ohms)
        s11_db[i] = 20 * np.log10(abs(gamma_in)) if abs(gamma_in) > 0 else -100
    
    return {
        'freq_ghz': freq_range_ghz,
        's11_db': s11_db,
        'freq_projeto_ghz': freq_projeto_hz / 1e9
    }

def calcular_casamento_lc(frequencia_hz, z0_ohms, zl_ohms):
    """
    Calcula os valores de Indutor (L) e Capacitor (C) para uma rede de casamento L.
    Retorna uma lista de todas as soluções possíveis.
    """
    if np.isclose(zl_ohms, z0_ohms):
        return {'solucoes': None, 'mensagem': "A impedância da carga já está casada."}

    w = 2 * np.pi * frequencia_hz  # Frequência angular
    solucoes = []

    RL = zl_ohms.real
    XL = zl_ohms.imag

    # --- Caso 1: RL > Z0 (Topologias com elemento paralelo primeiro na carga) ---
    if RL > z0_ohms:
        # Solução 1.1: Shunt C, Série L
        B_p = (XL + np.sqrt(RL / z0_ohms) * np.sqrt(RL**2 + XL**2 - z0_ohms * RL)) / (RL**2 + XL**2)
        Z_int = 1/(1/zl_ohms + 1j*B_p)
        X_s_final = -Z_int.imag
        if B_p > 0 and X_s_final > 0:
            C_shunt = B_p / w
            L_serie = X_s_final / w
            solucoes.append({'tipo': 'Shunt C, Série L', 'C': C_shunt, 'L': L_serie, 'z_int': Z_int})

        # Solução 1.2: Shunt L, Série C
        B_p = (XL - np.sqrt(RL / z0_ohms) * np.sqrt(RL**2 + XL**2 - z0_ohms * RL)) / (RL**2 + XL**2)
        Z_int = 1/(1/zl_ohms + 1j*B_p)
        X_s_final = -Z_int.imag
        if B_p < 0 and X_s_final < 0:
            L_shunt = -1 / (w * B_p)
            C_serie = -1 / (w * X_s_final)
            solucoes.append({'tipo': 'Shunt L, Série C', 'L': L_shunt, 'C': C_serie, 'z_int': Z_int})
    
    # --- Caso 2: RL < Z0 (Topologias com elemento em série primeiro na carga) ---
    if RL < z0_ohms and RL > 0:
        # Solução 2.1: Série L, Shunt C
        X_s = np.sqrt(RL * (z0_ohms - RL)) - XL
        Z_int = zl_ohms + 1j*X_s
        Y_int = 1/Z_int
        B_p_final = -Y_int.imag
        if B_p_final > 0:
            L_serie = X_s / w if X_s >= 0 else None # Permitir L=0
            C_shunt = B_p_final / w
            if L_serie is not None:
                solucoes.append({'tipo': 'Série L, Shunt C', 'L': L_serie, 'C': C_shunt, 'z_int': Z_int})

        # Solução 2.2: Série C, Shunt L
        X_s = -np.sqrt(RL * (z0_ohms - RL)) - XL
        Z_int = zl_ohms + 1j*X_s
        Y_int = 1/Z_int
        B_p_final = -Y_int.imag
        if B_p_final < 0:
            C_serie = -1/(w*X_s) if X_s < 0 else None
            L_shunt = -1/(w*B_p_final)
            if C_serie is not None:
                solucoes.append({'tipo': 'Série C, Shunt L', 'L': L_shunt, 'C': C_serie, 'z_int': Z_int})

    if not solucoes:
        return {'solucoes': None, 'mensagem': "Não foi encontrada uma solução de rede L simples.\nVerifique se RL > 0."}

    return {'solucoes': solucoes, 'mensagem': f"{len(solucoes)} solução(ões) encontrada(s)."}


# =============================================================================
# CLASSE PARA A CARTA DE SMITH INTERATIVA (VERSÃO REVISADA)
# =============================================================================

class SmithChart:
    def __init__(self, fig, ax):
        self.fig = fig
        self.ax = ax
        self.setup_chart()
        
    def _draw_wavelength_scales(self):
        radius = 1.08
        self.ax.text(0, radius + 0.1, "WAVELENGTHS TOWARD GENERATOR",
                             ha='center', va='center', fontsize=9, color='cyan', alpha=0.7)

        for i in range(51):
            val = i * 0.01
            angle_deg = 180 - (val * 720)
            angle_rad = np.deg2rad(angle_deg)
            
            is_major_tick = (i % 5 == 0)
            tick_len = 0.03 if is_major_tick else 0.015
            
            x_start, y_start = 1.0, 1.0
            x_end, y_end = (1 + tick_len), (1 + tick_len)

            self.ax.plot([x_start * np.cos(angle_rad), x_end * np.cos(angle_rad)],
                         [y_start * np.sin(angle_rad), y_end * np.sin(angle_rad)],
                         color='white', linewidth=1)
            
            if is_major_tick:
                label_rad = 1.05 + tick_len
                lx, ly = label_rad * np.cos(angle_rad), label_rad * np.sin(angle_rad)
                rot = angle_deg + 90 if angle_deg < 0 else angle_deg - 90
                self.ax.text(lx, ly, f'{val:.2f}', ha='center', va='center',
                             rotation=rot, fontsize=7)

    def setup_chart(self):
        self.ax.clear()
        self.ax.set_xlim(-1.4, 1.4)
        self.ax.set_ylim(-1.4, 1.4)
        self.ax.set_aspect('equal')
        self.ax.set_title('Carta de Smith Profissional', fontsize=12, fontweight='bold')
        
        outer_circle = Circle((0, 0), 1, fill=False, color='white', linewidth=2.0, zorder=10)
        self.ax.add_patch(outer_circle)
        
        grid_color = '#404040'
        r_values = [0.2, 0.5, 1.0, 2.0, 5.0]
        for r in r_values:
            center_x = r / (1 + r)
            radius = 1 / (1 + r)
            circle = Circle((center_x, 0), radius, fill=False, color=grid_color, alpha=0.9, linewidth=1.5)
            self.ax.add_patch(circle)
        
        x_values = [0.2, 0.5, 1.0, 2.0, 5.0, -0.2, -0.5, -1.0, -2.0, -5.0]
        for x in x_values:
            if x != 0:
                center_y = 1 / x
                radius = abs(1 / x)
                arc = patches.Arc((1, center_y), 2*radius, 2*radius, theta1=180 if x>0 else 0, theta2=0 if x>0 else 180, 
                                  fill=False, color=grid_color, alpha=0.9, linewidth=1.5)
                self.ax.add_patch(arc)

        self.ax.axhline(y=0, color=grid_color, linewidth=1.5, alpha=0.9)
        self.ax.axis('off')
        
        self._draw_wavelength_scales()

    def z_to_gamma(self, z_norm):
        return (z_norm - 1) / (z_norm + 1)
    
    def gamma_to_z(self, gamma):
        return (1 + gamma) / (1 - gamma)
    
    def y_to_gamma(self, y_norm):
        return (1 - y_norm) / (1 + y_norm)
    
    def rotate_gamma(self, gamma, beta_l):
        return gamma * np.exp(-2j * beta_l)
    
    def plot_impedance_matching(self, zl_norm, zl_ohms, d, l, lambda_0):
        self.setup_chart()
        beta = 2 * np.pi / lambda_0
        
        gamma_l = self.z_to_gamma(zl_norm)
        label_texto = f'ZL = {zl_ohms.real:.2f} + {zl_ohms.imag:+.2f}j Ω\n(zl = {zl_norm:.2f})'
        self.ax.plot(gamma_l.real, gamma_l.imag, 'o', c='#FF4136', markersize=10, label=label_texto, zorder=15)
        
        vswr = (1 + abs(gamma_l)) / (1 - abs(gamma_l)) if abs(gamma_l) < 1 else float('inf')
        if vswr != float('inf'):
            vswr_circle = Circle((0, 0), abs(gamma_l), fill=False, color='#FF851B', linewidth=2.5, alpha=0.9, linestyle='--', zorder=12)
            self.ax.add_patch(vswr_circle)
        
        gamma_d = self.rotate_gamma(gamma_l, beta * d)
        start_angle = np.angle(gamma_l, deg=True)
        end_angle = np.angle(gamma_d, deg=True)
        arc_path = patches.Arc((0, 0), 2*abs(gamma_l), 2*abs(gamma_l), angle=0,
                               theta1=end_angle, theta2=start_angle,
                               color='#2ECC40', linewidth=4, alpha=0.9,
                               label=f'Linha d ({d/lambda_0:.3f}λ)', zorder=13)
        self.ax.add_patch(arc_path)
        
        y_d = 1 / self.gamma_to_z(gamma_d)
        gamma_y_d = -gamma_d
        self.ax.plot(gamma_y_d.real, gamma_y_d.imag, 'D', c='#0074D9', markersize=9, label=f'Yd = {y_d:.2f}', zorder=15)
        
        if abs(y_d.imag) > 0.01:
            self.ax.plot([gamma_y_d.real, 0], [gamma_y_d.imag, 0], 
                         '#B10DC9', linewidth=4, alpha=0.9, marker='>', markersize=8,
                         label=f'Stub l ({l/lambda_0:.3f}λ)', zorder=14)
        
        self.ax.plot(0, 0, '*', c='yellow', markersize=18, label='Casamento (Z0)', zorder=16)
        
        info_text = f"λ = {lambda_0*100:.2f} cm\n"
        info_text += f"d = {d/lambda_0:.3f}λ ({d*100:.2f} cm)\n"
        info_text += f"l = {l/lambda_0:.3f}λ ({l*100:.2f} cm)\n"
        if vswr != float('inf'): info_text += f"VSWR = {vswr:.2f}\n"
        info_text += f"|Γ| = {abs(gamma_l):.3f}"
        
        self.ax.text(1.05, 0.5, info_text, fontsize=9, transform=self.ax.transAxes,
                     bbox=dict(boxstyle="round,pad=0.5", fc="#333", alpha=0.9, ec='cyan'),
                     verticalalignment='center', horizontalalignment='left')
        
        self.ax.legend(loc='upper left', bbox_to_anchor=(-0.3, 1.1), fontsize=8, framealpha=0.7)

    def plot_lc_matching(self, zl_norm, z_int_norm, zl_ohms, solucao):
        self.setup_chart()
        
        gamma_l = self.z_to_gamma(zl_norm)
        label_texto = f'ZL = {zl_ohms.real:.2f} + {zl_ohms.imag:+.2f}j Ω\n(zl = {zl_norm:.2f})'
        self.ax.plot(gamma_l.real, gamma_l.imag, 'o', c='#FF4136', markersize=10, label=label_texto, zorder=15)

        gamma_int = self.z_to_gamma(z_int_norm)
        self.ax.plot(gamma_int.real, gamma_int.imag, 'D', c='#0074D9', markersize=9, label=f'Ponto Intermediário\n(z_int = {z_int_norm:.2f})', zorder=15)

        self.ax.plot(0, 0, '*', c='yellow', markersize=18, label='Casamento (Z0)', zorder=16)

        tipo = solucao['tipo']
        # Desenha arcos de componentes
        if 'Série' in tipo.split(',')[0]:
            r_const = zl_norm.real
            center_x = r_const / (r_const + 1)
            radius = 1 / (r_const + 1)
            gamma_center = complex(center_x, 0)
            theta1 = np.angle(gamma_l - gamma_center, deg=True)
            theta2 = np.angle(gamma_int - gamma_center, deg=True)
            arc1 = patches.Arc((center_x, 0), 2*radius, 2*radius, angle=0, theta1=theta1, theta2=theta2, color='#2ECC40', linewidth=4, alpha=0.9, zorder=13)
            self.ax.add_patch(arc1)
        else: # Shunt
            y_norm = 1 / zl_norm
            g_const = y_norm.real
            center_x = -g_const / (g_const + 1)
            radius = 1 / (g_const + 1)
            gamma_center = complex(center_x, 0)
            theta1 = np.angle(-gamma_l - gamma_center, deg=True)
            theta2 = np.angle(-gamma_int - gamma_center, deg=True)
            arc1 = patches.Arc((center_x, 0), 2*radius, 2*radius, angle=0, theta1=theta1, theta2=theta2, color='#2ECC40', linewidth=4, alpha=0.9, zorder=13)
            self.ax.add_patch(arc1)
        
        if 'Série' in tipo.split(',')[1]:
            r_const = z_int_norm.real
            center_x = r_const / (r_const + 1)
            radius = 1 / (r_const + 1)
            gamma_center = complex(center_x, 0)
            theta1 = np.angle(gamma_int - gamma_center, deg=True)
            theta2 = np.angle(0 - gamma_center, deg=True)
            arc2 = patches.Arc((center_x, 0), 2*radius, 2*radius, angle=0, theta1=theta1, theta2=theta2, color='#B10DC9', linewidth=4, alpha=0.9, zorder=14)
            self.ax.add_patch(arc2)
        else: # Shunt
            y_int_norm = 1 / z_int_norm
            g_const = y_int_norm.real
            center_x = -g_const / (g_const + 1)
            radius = 1 / (g_const + 1)
            gamma_center = complex(center_x, 0)
            theta1 = np.angle(-gamma_int - gamma_center, deg=True)
            theta2 = np.angle(0 - gamma_center, deg=True)
            arc2 = patches.Arc((center_x, 0), 2*radius, 2*radius, angle=0, theta1=theta1, theta2=theta2, color='#B10DC9', linewidth=4, alpha=0.9, zorder=14)
            self.ax.add_patch(arc2)
        
        info_text = f"Topologia: {tipo}\n"
        if 'L' in solucao: info_text += f"Indutor: {solucao['L']*1e9:.2f} nH\n"
        if 'C' in solucao: info_text += f"Capacitor: {solucao['C']*1e12:.2f} pF"
        
        self.ax.text(1.05, 0.5, info_text, fontsize=9, transform=self.ax.transAxes,
                     bbox=dict(boxstyle="round,pad=0.5", fc="#333", alpha=0.9, ec='cyan'),
                     verticalalignment='center', horizontalalignment='left')

        self.ax.legend(loc='upper left', bbox_to_anchor=(-0.3, 1.1), fontsize=8, framealpha=0.7)

    def clear_chart(self):
        self.setup_chart()
        
# =============================================================================
# CLASSE DA APLICAÇÃO TKINTER (LÓGICA DO FRONT-END)
# =============================================================================
class CalculatorApp:
    def __init__(self, master):
        self.master = master
        master.title("Calculadora de Casamento de Impedâncias")
        master.resizable(True, True)

        self.stub_results = None
        self.lc_results = None
        self.microstrip_dims = {}

        self.notebook = tb.Notebook(master)
        self.notebook.pack(fill="both", expand=True, padx=10, pady=10)
        
        main_tab = tb.Frame(self.notebook, padding=10)
        self.notebook.add(main_tab, text="Casamento com Stub")
        
        lc_tab = tb.Frame(self.notebook, padding=10)
        self.notebook.add(lc_tab, text="Casamento LC")
        
        bw_tab = tb.Frame(self.notebook, padding=10)
        self.notebook.add(bw_tab, text="Análise de Largura de Banda")
        
        visual_tab = tb.Frame(self.notebook, padding=10)
        self.notebook.add(visual_tab, text="Visualização da Microfita")

        self._setup_main_tab(main_tab)
        self._setup_lc_tab(lc_tab)
        self._setup_bw_tab(bw_tab)
        self._setup_visualization_tab(visual_tab)
    
    def _setup_main_tab(self, main_frame):
        main_frame.grid_rowconfigure(1, weight=1)
        main_frame.grid_columnconfigure(1, weight=1)

        input_frame = tb.Labelframe(main_frame, text="Parâmetros de Entrada", padding=10)
        input_frame.grid(row=0, column=0, columnspan=2, sticky="ew")

        results_frame = tb.Labelframe(main_frame, text="Resultados", padding=10)
        results_frame.grid(row=1, column=0, sticky="nsew", pady=5)
        
        chart_frame = tb.Labelframe(main_frame, text="Carta de Smith Interativa", padding=10)
        chart_frame.grid(row=1, column=1, sticky="nsew", padx=(5, 0), pady=5)
        
        self.entries = {}
        input_fields = {
            "Cálculo do Stub": {
                "Frequência (GHz):": "freq",
                "Impedância Z0 (Ω):": "z0",
                "Carga ZL (Real) (Ω):": "zl_real",
                "Carga ZL (Imag) (Ω):": "zl_imag"
            },
            "Dimensionamento da Microfita (Opcional)": {
                "Altura do Substrato h (mm):": "h",
                "Constante Dielétrica εr:": "er"
            }
        }
        
        row_counter = 0
        for section_title, fields in input_fields.items():
            tb.Label(input_frame, text=section_title, font=("Helvetica", 10, "bold")).grid(row=row_counter, column=0, columnspan=2, sticky="w", pady=(10, 2))
            row_counter += 1
            for label_text, key in fields.items():
                label = tb.Label(input_frame, text=label_text)
                label.grid(row=row_counter, column=0, sticky="w", padx=5, pady=2)
                entry = tb.Entry(input_frame, width=20)
                entry.grid(row=row_counter, column=1, sticky="ew", padx=5, pady=2)
                self.entries[key] = entry
                row_counter += 1

        button_frame = tb.Frame(main_frame)
        button_frame.grid(row=2, column=0, columnspan=2, pady=10)
        
        calc_button = tb.Button(button_frame, text="Calcular Stub", command=self.calculate)
        calc_button.grid(row=0, column=0, padx=5)

        clear_button = tb.Button(button_frame, text="Limpar", command=self.clear_fields)
        clear_button.grid(row=0, column=1, padx=5)
        
        demo_button = tb.Button(button_frame, text="Demo", command=self.load_demo)
        demo_button.grid(row=0, column=2, padx=5)

        text_bg = self.master.cget("background")
        self.results_text = tk.Text(results_frame, height=12, width=50, wrap="word", state="disabled", font=("Courier", 9), bg=text_bg, fg="white", highlightthickness=0, bd=0)
        self.results_text.grid(row=0, column=0, sticky="nsew")
        scrollbar = tb.Scrollbar(results_frame, orient="vertical", command=self.results_text.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.results_text['yscrollcommand'] = scrollbar.set
        
        results_frame.grid_rowconfigure(0, weight=1)
        results_frame.grid_columnconfigure(0, weight=1)
        
        plt.style.use('dark_background')
        self.fig, self.ax = plt.subplots(figsize=(7, 7), facecolor='#2b2b2b', constrained_layout=True)
        self.ax.set_facecolor('#2b2b2b')
        
        self.smith_chart = SmithChart(self.fig, self.ax)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=chart_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        
        self.status_label = tb.Label(chart_frame, text="Mova o mouse sobre a carta para ver valores", font=("Courier", 8))
        self.status_label.grid(row=1, column=0, pady=5)
        
        chart_frame.grid_rowconfigure(0, weight=1)
        chart_frame.grid_columnconfigure(0, weight=1)

    def _setup_lc_tab(self, lc_frame):
        # <<< MODIFIED: Changed layout to support diagram
        lc_frame.grid_columnconfigure(0, weight=1)
        lc_frame.grid_columnconfigure(1, weight=2) # Give more space to diagram
        lc_frame.grid_rowconfigure(1, weight=1)

        input_frame = tb.Labelframe(lc_frame, text="Parâmetros de Entrada", padding=10)
        # Place it across both columns at the top
        input_frame.grid(row=0, column=0, columnspan=2, sticky="ew", padx=5, pady=5)

        results_frame = tb.Labelframe(lc_frame, text="Soluções Encontradas", padding=10)
        # Place it in the bottom-left cell
        results_frame.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)
        results_frame.grid_rowconfigure(0, weight=1)
        results_frame.grid_columnconfigure(0, weight=1)
        
        # <<< NEW: Frame for the LC schematic diagram
        diagram_frame = tb.Labelframe(lc_frame, text="Diagrama da Rede", padding=10)
        diagram_frame.grid(row=1, column=1, sticky="nsew", padx=5, pady=5)

        self.lc_entries = {}
        input_fields = {
            "Frequência (GHz):": "freq", "Impedância Z0 (Ω):": "z0",
            "Carga ZL (Real) (Ω):": "zl_real", "Carga ZL (Imag) (Ω):": "zl_imag"
        }
        
        for i, (label_text, key) in enumerate(input_fields.items()):
            tb.Label(input_frame, text=label_text).grid(row=i, column=0, sticky="w", padx=5, pady=3)
            entry = tb.Entry(input_frame, width=20)
            entry.grid(row=i, column=1, sticky="ew", padx=5, pady=3)
            self.lc_entries[key] = entry
            input_frame.grid_columnconfigure(1, weight=1)
        
        calc_button = tb.Button(input_frame, text="Calcular Redes LC", command=self.calculate_lc)
        calc_button.grid(row=0, column=2, rowspan=len(input_fields), padx=20)
        
        self.lc_results_text = tk.Text(results_frame, height=10, width=50, wrap="word", state="disabled", font=("Courier", 10))
        self.lc_results_text.grid(row=0, column=0, sticky="nsew")
        
        scrollbar = tb.Scrollbar(results_frame, orient="vertical", command=self.lc_results_text.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.lc_results_text['yscrollcommand'] = scrollbar.set
        
        self.plot_buttons_frame = tb.Frame(results_frame)
        self.plot_buttons_frame.grid(row=1, column=0, columnspan=2, pady=5)
        
        # <<< NEW: Setup Matplotlib canvas for the LC schematic
        self.lc_fig, self.lc_ax = plt.subplots(figsize=(6, 4), facecolor='#2b2b2b', constrained_layout=True)
        self.lc_canvas = FigureCanvasTkAgg(self.lc_fig, master=diagram_frame)
        self.lc_canvas.get_tk_widget().pack(fill="both", expand=True)
        self._clear_lc_diagram()


    def _setup_bw_tab(self, bw_tab):
        bw_frame = tb.Frame(bw_tab)
        bw_frame.pack(fill="both", expand=True)
        
        bw_control_frame = tb.Labelframe(bw_frame, text="Parâmetros de Varredura", padding=10)
        bw_control_frame.pack(fill="x", padx=10, pady=10)
        
        self.freq_entries = {}
        freq_fields = {
            "Frequência Inicial (GHz):": "freq_min",
            "Frequência Final (GHz):": "freq_max",
            "Passo (GHz):": "freq_step"
        }
        
        for i, (label_text, key) in enumerate(freq_fields.items()):
            label = tb.Label(bw_control_frame, text=label_text)
            label.grid(row=0, column=i*2, sticky="w", padx=5, pady=2)
            entry = tb.Entry(bw_control_frame, width=10)
            entry.grid(row=0, column=i*2+1, sticky="ew", padx=5, pady=2)
            self.freq_entries[key] = entry
        
        self.s11_button = tb.Button(bw_control_frame, text="Gerar Gráfico S11", command=self.plot_s11)
        self.s11_button.grid(row=0, column=len(freq_fields)*2, padx=15)
        
        s11_chart_frame = tb.Labelframe(bw_frame, text="Gráfico de Perda de Retorno (S11) vs. Frequência", padding=10)
        s11_chart_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        self.s11_fig, self.s11_ax = plt.subplots(figsize=(8, 5), facecolor='#2b2b2b', constrained_layout=True)
        self.s11_ax.set_facecolor('#2b2b2b')
        self.s11_ax.set_xlabel('Frequência (GHz)')
        self.s11_ax.set_ylabel('S11 (dB)')
        self.s11_ax.grid(True, alpha=0.3)
        
        self.s11_canvas = FigureCanvasTkAgg(self.s11_fig, master=s11_chart_frame)
        self.s11_canvas.draw()
        self.s11_canvas.get_tk_widget().pack(fill="both", expand=True)
        
        self.bw_info_label = tb.Label(s11_chart_frame, text="A largura de banda será exibida aqui", font=("Courier", 9))
        self.bw_info_label.pack(pady=5)
        
        self.freq_entries["freq_min"].insert(0, "2.0")
        self.freq_entries["freq_max"].insert(0, "3.0")
        self.freq_entries["freq_step"].insert(0, "0.01")

    def _setup_visualization_tab(self, visual_tab):
        self.vis_fig, (self.vis_ax1, self.vis_ax2) = plt.subplots(
            2, 1, figsize=(8, 8), facecolor='#2b2b2b', 
            gridspec_kw={'height_ratios': [1, 2], 'hspace': 0.4},
            constrained_layout=True
        )
        self.vis_fig.set_facecolor('#2b2b2b')
        
        self.vis_canvas = FigureCanvasTkAgg(self.vis_fig, master=visual_tab)
        self.vis_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        self._clear_visualization_plot()

    def _clear_visualization_plot(self, message="Calcule uma solução para ver a visualização da microfita."):
        for ax in [self.vis_ax1, self.vis_ax2]:
            ax.clear()
            ax.set_facecolor('#2b2b2b')
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)

        self.vis_ax1.set_title("Visão de Corte Transversal", fontweight='bold')
        self.vis_ax2.set_title("Visão Superior (Layout)", fontweight='bold')

        self.vis_ax2.text(0.5, 0.5, message,
                          ha='center', va='center', fontsize=12, color='gray', transform=self.vis_ax2.transAxes)
        self.vis_canvas.draw()
    
    def _clear_lc_diagram(self, message="Selecione uma solução para ver o diagrama."):
        self.lc_ax.clear()
        self.lc_ax.set_facecolor('#2b2b2b')
        self.lc_ax.set_xlim(0, 10.5)
        self.lc_ax.set_ylim(0, 6)
        for spine in self.lc_ax.spines.values():
            spine.set_visible(False)
        self.lc_ax.set_xticks([])
        self.lc_ax.set_yticks([])
        self.lc_ax.text(5.25, 3, message, ha='center', va='center', fontsize=12, color='gray')
        self.lc_canvas.draw()

    # =============================================================================
    # >>> INÍCIO DA FUNÇÃO CORRIGIDA <<<
    # =============================================================================
    def draw_lc_schematic(self, solution):
        self._clear_lc_diagram(message="") # Limpa o diagrama anterior
        ax = self.lc_ax
        tipo = solution['tipo']
        ax.set_title(f"Topologia: {tipo}", fontweight='bold')
        
        # --- Parâmetros de Desenho ---
        line_color, font_size, comp_color, load_color = 'white', 10, 'cyan', 'lime'
        
        # --- Coordenadas ---
        y_main, y_gnd = 4.0, 1.0
        x_in_start, x_in_end = 1.0, 2.5
        x_source_comp, x_load_comp = 4.0, 7.0
        x_load_box_start = 8.5
        
        # --- Desenha Elementos Comuns ---
        # Linha de Terra
        ax.plot([x_in_start, x_load_box_start + 0.75], [y_gnd, y_gnd], color=line_color, lw=2)
        # Entrada (Zin)
        ax.text(x_in_start - 0.5, y_main, r'$Z_{in}$', ha='center', va='center', fontsize=font_size, color=load_color)
        ax.plot([x_in_start, x_in_end], [y_main, y_main], color=line_color, lw=2)
        
        # Carga (ZL) - CORRIGIDO para ser um quadrado com conexões explícitas
        load_side = 1.5
        load_box_y = (y_main + y_gnd) / 2 - load_side / 2
        load_center_x = x_load_box_start + load_side / 2
        ax.add_patch(patches.Rectangle((x_load_box_start, load_box_y), load_side, load_side,
                                        facecolor='none', edgecolor=load_color, lw=2))
        ax.text(load_center_x, (y_main + y_gnd) / 2, r'$Z_L$',
                ha='center', va='center', fontsize=font_size, color=load_color)
        # Conexões verticais para a carga
        ax.plot([load_center_x, load_center_x], [y_main, load_box_y + load_side], color=line_color, lw=2)
        ax.plot([load_center_x, load_center_x], [load_box_y, y_gnd], color=line_color, lw=2)

        # A convenção da topologia é: 'Componente Perto da Carga, Componente Perto da Fonte'
        if tipo == 'Shunt C, Série L':
            # --- Componentes ---
            # Shunt C (perto da carga, no ponto da linha principal)
            w, gap = 0.8, 0.2
            y_top, y_bot = y_main - 0.3, y_main - 0.3 - gap
            ax.plot([x_load_comp, x_load_comp], [y_main, y_top], color=line_color, lw=2)
            ax.plot([x_load_comp - w/2, x_load_comp + w/2], [y_top, y_top], color=line_color, lw=2)
            ax.plot([x_load_comp - w/2, x_load_comp + w/2], [y_bot, y_bot], color=line_color, lw=2)
            ax.plot([x_load_comp, x_load_comp], [y_bot, y_gnd], color=line_color, lw=2)
            ax.text(x_load_comp, y_gnd - 0.5, f"C={solution['C']*1e12:.2f} pF", ha='center', va='center', fontsize=font_size, color=comp_color)
            
            # Série L (perto da fonte, quebra a linha) - VISUAL MELHORADO
            r, sep, n = 0.2, 0.3, 4
            width = (n - 1) * sep
            half_w = (width + 2*r) / 2
            start_x = x_source_comp - width/2
            for i in range(n):
                ax.add_patch(patches.Arc((start_x + i*sep, y_main), 2*r, 2*r, theta1=180, theta2=0, ec=line_color, lw=2))
            ax.text(x_source_comp, y_main + 1.2, f"L={solution['L']*1e9:.2f} nH", ha='center', va='center', fontsize=font_size, color=comp_color)
            
            # --- Conexões ---
            ax.plot([x_in_end, x_source_comp - half_w], [y_main, y_main], color=line_color, lw=2)
            ax.plot([x_source_comp + half_w, load_center_x], [y_main, y_main], color=line_color, lw=2)
            
        elif tipo == 'Shunt L, Série C':
            # --- Componentes ---
            # Shunt L (perto da carga)
            r, n = 0.25, 3
            y_start = y_main - 0.3; height = n * r * 2
            ax.plot([x_load_comp, x_load_comp], [y_main, y_start], color=line_color, lw=2)
            for i in range(n):
                ax.add_patch(patches.Arc((x_load_comp, y_start - r - i*2*r), 2*r, 2*r, theta1=270, theta2=90, ec=line_color, lw=2))
            ax.plot([x_load_comp, x_load_comp], [y_start - height, y_gnd], color=line_color, lw=2)
            ax.text(x_load_comp, y_gnd - 0.5, f"L={solution['L']*1e9:.2f} nH", ha='center', va='center', fontsize=font_size, color=comp_color)
            
            # Série C (perto da fonte)
            gap, h, half_w = 0.2, 1.0, 0.2
            ax.plot([x_source_comp-gap, x_source_comp-gap], [y_main-h/2, y_main+h/2], color=line_color, lw=2)
            ax.plot([x_source_comp+gap, x_source_comp+gap], [y_main-h/2, y_main+h/2], color=line_color, lw=2)
            ax.text(x_source_comp, y_main+1.2, f"C={solution['C']*1e12:.2f} pF", ha='center', va='center', fontsize=font_size, color=comp_color)
            
            # --- Conexões ---
            ax.plot([x_in_end, x_source_comp - half_w], [y_main, y_main], color=line_color, lw=2)
            ax.plot([x_source_comp + half_w, load_center_x], [y_main, y_main], color=line_color, lw=2)
        
        elif tipo == 'Série L, Shunt C':
            # --- Componentes ---
            # Série L (perto da carga) - VISUAL MELHORADO
            r, sep, n = 0.2, 0.3, 4
            width = (n - 1) * sep
            half_w = (width + 2*r) / 2
            start_x = x_load_comp - width/2
            for i in range(n):
                ax.add_patch(patches.Arc((start_x + i*sep, y_main), 2*r, 2*r, theta1=180, theta2=0, ec=line_color, lw=2))
            ax.text(x_load_comp, y_main+1.2, f"L={solution['L']*1e9:.2f} nH", ha='center', va='center', fontsize=font_size, color=comp_color)

            # Shunt C (perto da fonte)
            w, gap = 0.8, 0.2
            y_top, y_bot = y_main - 0.3, y_main - 0.3 - gap
            ax.plot([x_source_comp, x_source_comp], [y_main, y_top], color=line_color, lw=2)
            ax.plot([x_source_comp - w/2, x_source_comp + w/2], [y_top, y_top], color=line_color, lw=2)
            ax.plot([x_source_comp - w/2, x_source_comp + w/2], [y_bot, y_bot], color=line_color, lw=2)
            ax.plot([x_source_comp, x_source_comp], [y_bot, y_gnd], color=line_color, lw=2)
            ax.text(x_source_comp, y_gnd - 0.5, f"C={solution['C']*1e12:.2f} pF", ha='center', va='center', fontsize=font_size, color=comp_color)
            
            # --- Conexões ---
            ax.plot([x_in_end, x_load_comp - half_w], [y_main, y_main], color=line_color, lw=2)
            ax.plot([x_load_comp + half_w, load_center_x], [y_main, y_main], color=line_color, lw=2)

        elif tipo == 'Série C, Shunt L':
            # --- Componentes ---
            # Série C (perto da carga)
            gap, h, half_w = 0.2, 1.0, 0.2
            ax.plot([x_load_comp - gap, x_load_comp - gap], [y_main - h/2, y_main + h/2], color=line_color, lw=2)
            ax.plot([x_load_comp + gap, x_load_comp + gap], [y_main - h/2, y_main + h/2], color=line_color, lw=2)
            ax.text(x_load_comp, y_main+1.2, f"C={solution['C']*1e12:.2f} pF", ha='center', va='center', fontsize=font_size, color=comp_color)
            
            # Shunt L (perto da fonte)
            r, n = 0.25, 3
            y_start = y_main - 0.3; height = n * r * 2
            ax.plot([x_source_comp, x_source_comp], [y_main, y_start], color=line_color, lw=2)
            for i in range(n):
                ax.add_patch(patches.Arc((x_source_comp, y_start - r - i*2*r), 2*r, 2*r, theta1=270, theta2=90, ec=line_color, lw=2))
            ax.plot([x_source_comp, x_source_comp], [y_start - height, y_gnd], color=line_color, lw=2)
            ax.text(x_source_comp, y_gnd - 0.5, f"L={solution['L']*1e9:.2f} nH", ha='center', va='center', fontsize=font_size, color=comp_color)
            
            # --- Conexões ---
            ax.plot([x_in_end, x_load_comp - half_w], [y_main, y_main], color=line_color, lw=2)
            ax.plot([x_load_comp + half_w, load_center_x], [y_main, y_main], color=line_color, lw=2)

        self.lc_canvas.draw()
    # =============================================================================
    # >>> FIM DA FUNÇÃO CORRIGIDA <<<
    # =============================================================================

    def draw_microstrip_visualization(self):
        self._clear_visualization_plot(message="")
        
        w_mm = self.microstrip_dims.get('w_mm')
        h_mm = self.microstrip_dims.get('h_mm')
        er = self.microstrip_dims.get('er')
        sol = self.stub_results['solucoes']
        d_m, l_m = sol['d'], sol['l']

        if w_mm and h_mm and er:
            self.vis_ax1.set_title("Visão de Corte Transversal", fontweight='bold')
            ground_h, strip_h, total_w = h_mm * 0.2, h_mm * 0.15, w_mm * 5
            ground = patches.Rectangle((-total_w/2, -ground_h), total_w, ground_h, color='gray')
            substrato = patches.Rectangle((-total_w/2, 0), total_w, h_mm, color='darkgreen', alpha=0.8)
            strip = patches.Rectangle((-w_mm/2, h_mm), w_mm, strip_h, color='orange')
            self.vis_ax1.add_patch(ground); self.vis_ax1.add_patch(substrato); self.vis_ax1.add_patch(strip)
            self.vis_ax1.text(0, h_mm/2, f'εr = {er}', color='white', ha='center', va='center', fontsize=10, fontweight='bold')
            self.vis_ax1.annotate('', xy=(total_w/2 * 0.8, 0), xytext=(total_w/2 * 0.8, h_mm), arrowprops=dict(arrowstyle='<->', color='cyan'))
            self.vis_ax1.text(total_w/2 * 0.85, h_mm/2, f'h = {h_mm:.2f} mm', color='cyan', ha='left', va='center', rotation=90)
            self.vis_ax1.annotate('', xy=(-w_mm/2, h_mm + strip_h*1.5), xytext=(w_mm/2, h_mm + strip_h*1.5), arrowprops=dict(arrowstyle='<->', color='yellow'))
            self.vis_ax1.text(0, h_mm + strip_h*2, f'W = {w_mm:.3f} mm', color='yellow', ha='center', va='bottom')
            self.vis_ax1.set_xlim(-total_w/1.8, total_w/1.8); self.vis_ax1.set_ylim(-ground_h, h_mm + strip_h * 4)
            self.vis_ax1.set_aspect('equal', adjustable='box')
        else:
            self.vis_ax1.text(0.5, 0.5, "Forneça 'h' e 'εr' para ver o corte transversal.",
                              ha='center', va='center', fontsize=10, color='gray', transform=self.vis_ax1.transAxes)

        self.vis_ax2.set_title("Visão Superior (Layout)", fontweight='bold')
        w_vis = w_mm/1000 if w_mm else max(d_m, l_m) * 0.1
        load_pos_x, stub_pos_x = 0, -d_m
        line_len = d_m + max(d_m, l_m) * 0.3
        main_line = patches.Rectangle((stub_pos_x - line_len + d_m, -w_vis/2), line_len, w_vis, color='orange')
        stub_line = patches.Rectangle((stub_pos_x - w_vis/2, -l_m), w_vis, l_m, color='orange')
        short = patches.Rectangle((stub_pos_x - w_vis, -l_m - w_vis*0.1), w_vis*2, w_vis*0.2, color='black')
        self.vis_ax2.add_patch(main_line); self.vis_ax2.add_patch(stub_line); self.vis_ax2.add_patch(short)
        self.vis_ax2.text(load_pos_x + w_vis, 0, 'Carga\n(ZL)', color='red', ha='left', va='center')
        self.vis_ax2.text(stub_pos_x - line_len + d_m - w_vis, 0, 'Fonte', color='lime', ha='right', va='center')
        self.vis_ax2.annotate('', xy=(stub_pos_x, w_vis*2), xytext=(load_pos_x, w_vis*2), arrowprops=dict(arrowstyle='<->', color='cyan'))
        self.vis_ax2.text((stub_pos_x + load_pos_x)/2, w_vis*2.5, f'd = {d_m*100:.2f} cm', color='cyan', ha='center', va='bottom')
        self.vis_ax2.annotate('', xy=(stub_pos_x + w_vis, -l_m), xytext=(stub_pos_x + w_vis, 0), arrowprops=dict(arrowstyle='<->', color='yellow'))
        self.vis_ax2.text(stub_pos_x + w_vis*1.5, -l_m/2, f'l = {l_m*100:.2f} cm', color='yellow', ha='left', va='center')
        self.vis_ax2.set_xlim(stub_pos_x - line_len + d_m - w_vis*2, load_pos_x + w_vis*4)
        self.vis_ax2.set_ylim(-l_m - w_vis*3, w_vis*6)
        self.vis_ax2.set_aspect('equal', adjustable='box')
        self.vis_canvas.draw()

    def on_mouse_move(self, event):
        if event.inaxes == self.ax and event.xdata is not None and event.ydata is not None:
            gamma = complex(event.xdata, event.ydata)
            if abs(gamma) <= 1.0:
                z_norm = self.smith_chart.gamma_to_z(gamma)
                y_norm = 1 / z_norm
                status_text = f"Γ = {gamma:.3f} | Z = {z_norm:.3f} | Y = {y_norm:.3f}"
                self.status_label.config(text=status_text)
            else:
                self.status_label.config(text="Fora da carta de Smith")
        else:
            self.status_label.config(text="Mova o mouse sobre a carta para ver valores")

    def clear_fields(self):
        for entry in self.entries.values(): entry.delete(0, tk.END)
        self.results_text.config(state="normal"); self.results_text.delete("1.0", tk.END); self.results_text.config(state="disabled")
        self.smith_chart.clear_chart(); self.canvas.draw()
        self.s11_ax.clear(); self.s11_ax.grid(True, alpha=0.3)
        self.s11_ax.set_xlabel('Frequência (GHz)'); self.s11_ax.set_ylabel('S11 (dB)')
        self.s11_canvas.draw(); self.bw_info_label.config(text="A largura de banda será exibida aqui")
        self._clear_visualization_plot()
        self._clear_lc_diagram()
        self.stub_results = None; self.microstrip_dims = {}

    def calculate(self):
        try:
            f_ghz = float(self.entries["freq"].get()); z0 = float(self.entries["z0"].get())
            zl_real = float(self.entries["zl_real"].get()); zl_imag = float(self.entries["zl_imag"].get())
            zl = complex(zl_real, zl_imag)
            h_str = self.entries["h"].get(); er_str = self.entries["er"].get()
            
            self.results_text.config(state="normal"); self.results_text.delete("1.0", tk.END)
            output_str = ""
            self.stub_results = calcular_stub_paralelo_curto(f_ghz * 1e9, z0, zl)
            
            if not self.stub_results['solucoes']:
                output_str += f"STUB: {self.stub_results['mensagem']}\n"
                self.smith_chart.clear_chart(); self.canvas.draw()
                self._clear_visualization_plot(message="Cálculo do stub falhou.\nVisualização indisponível.")
            else:
                sol = self.stub_results['solucoes']
                output_str += "--- RESULTADO DO STUB ---\n"
                output_str += f"  Distância (d): {sol['d']*100:.3f} cm\n"
                output_str += f"  Comprimento (l): {sol['l']*100:.3f} cm\n"
                output_str += f"  (d = {sol['d']/self.stub_results['lambda']:.3f} λ, l = {sol['l']/self.stub_results['lambda']:.3f} λ)\n"

            output_str += "-"*45 + "\n"
            
            self.microstrip_dims = {}
            if h_str and er_str:
                h_mm = float(h_str); er = float(er_str)
                largura_w_m = dimensionar_microstrip(z0, h_mm / 1000.0, er)
                if largura_w_m is not None:
                    output_str += "--- DIMENSIONAMENTO DA MICROFITA ---\n"
                    output_str += f"  Para Z0 = {z0:.1f} Ω, h = {h_mm:.3f} mm, εr = {er:.2f}\n"
                    output_str += f"  Largura da Linha (W) = {largura_w_m * 1000:.4f} mm\n"
                    self.microstrip_dims = {'w_mm': largura_w_m * 1000, 'h_mm': h_mm, 'er': er}
                else: output_str += "MICROFITA: Não foi possível calcular a largura.\n"
            else: output_str += "MICROFITA: Preencha 'h' e 'εr' para calcular.\n"

            self.results_text.insert("1.0", output_str); self.results_text.config(state="disabled")
            
            if self.stub_results['solucoes']:
                zl_norm = zl / z0; sol = self.stub_results['solucoes']
                self.smith_chart.plot_impedance_matching(zl_norm, zl, sol['d'], sol['l'], self.stub_results['lambda'])
                self.canvas.draw(); self.draw_microstrip_visualization(); self.notebook.select(3)
                self.freq_entries["freq_min"].delete(0, tk.END); self.freq_entries["freq_min"].insert(0, f"{max(0.5, f_ghz - 0.5):.1f}")
                self.freq_entries["freq_max"].delete(0, tk.END); self.freq_entries["freq_max"].insert(0, f"{f_ghz + 0.5:.1f}")
            else:
                self.smith_chart.clear_chart(); self.canvas.draw()
        except ValueError: Messagebox.show_error("Erro de Entrada", "Verifique se todos os campos numéricos foram preenchidos corretamente.")
        except Exception as e: Messagebox.show_error("Erro Inesperado", f"Ocorreu um erro: {e}")

    def calculate_lc(self):
        try:
            f_ghz = float(self.lc_entries["freq"].get())
            z0 = float(self.lc_entries["z0"].get())
            zl_real = float(self.lc_entries["zl_real"].get())
            zl_imag = float(self.lc_entries["zl_imag"].get())
            zl = complex(zl_real, zl_imag)
            
            self.lc_results = calcular_casamento_lc(f_ghz * 1e9, z0, zl)

            self.lc_results_text.config(state="normal")
            self.lc_results_text.delete("1.0", tk.END)
            
            for widget in self.plot_buttons_frame.winfo_children():
                widget.destroy()

            if not self.lc_results['solucoes']:
                self.lc_results_text.insert("1.0", self.lc_results['mensagem'])
                self._clear_lc_diagram()
            else:
                output_str = f"--- {self.lc_results['mensagem']} ---\n\n"
                for i, sol in enumerate(self.lc_results['solucoes']):
                    output_str += f"SOLUÇÃO #{i+1}: {sol['tipo']}\n"
                    if 'L' in sol: output_str += f"  L = {sol['L'] * 1e9:.3f} nH\n"
                    if 'C' in sol: output_str += f"  C = {sol['C'] * 1e12:.3f} pF\n"
                    output_str += "-"*30 + "\n"
                    
                    plot_btn = tb.Button(self.plot_buttons_frame, text=f"Plotar & Ver Solução #{i+1}", 
                                         command=lambda s=sol: self.plot_selected_lc(s, z0, zl))
                    plot_btn.pack(side="left", padx=5)

                self.lc_results_text.insert("1.0", output_str)
                first_solution = self.lc_results['solucoes'][0]
                self.plot_selected_lc(first_solution, z0, zl)
                self.notebook.select(1) # Switch to LC tab to show the user

        except ValueError:
            Messagebox.show_error("Erro de Entrada", "Verifique os campos numéricos.")
        except Exception as e:
            Messagebox.show_error("Erro Inesperado", f"Ocorreu um erro: {e}")
        finally:
            self.lc_results_text.config(state="disabled")

    def plot_selected_lc(self, solucao, z0, zl):
        zl_norm = zl / z0
        z_int_norm = solucao['z_int'] / z0
        
        # Plot on Smith Chart
        self.smith_chart.plot_lc_matching(zl_norm, z_int_norm, zl, solucao)
        self.canvas.draw()
        
        # Draw the schematic diagram
        self.draw_lc_schematic(solucao)
        
        # Go to the Smith Chart tab to show the result immediately
        self.notebook.select(0)

    def plot_s11(self):
        try:
            if not self.stub_results or not self.stub_results['solucoes']:
                Messagebox.show_warning("Aviso", "Primeiro calcule uma solução de stub na aba principal."); return
            
            freq_min = float(self.freq_entries["freq_min"].get()); freq_max = float(self.freq_entries["freq_max"].get()); freq_step = float(self.freq_entries["freq_step"].get())
            if freq_min >= freq_max or freq_step <= 0: Messagebox.show_error("Erro", "Verifique os valores da faixa de frequência."); return
            
            z0 = self.stub_results['z0_ohms']; zl = self.stub_results['zl_ohms']; solucao = self.stub_results['solucoes']
            solucao['frequencia_hz'] = self.stub_results['frequencia_hz']
            s11_results = calcular_s11_vs_frequencia(freq_min, freq_max, freq_step, z0, zl, solucao)
            
            if not s11_results: Messagebox.show_error("Erro", "Não foi possível calcular o S11."); return
            
            self.s11_ax.clear()
            self.s11_ax.plot(s11_results['freq_ghz'], s11_results['s11_db'], 'cyan', linewidth=2)
            self.s11_ax.axhline(-10, color='red', linestyle='--', alpha=0.7)
            self.s11_ax.text(freq_min, -9.5, '-10 dB', color='red', fontsize=8, ha='left')
            
            proj_freq = s11_results.get('freq_projeto_ghz', 0)
            if proj_freq > 0:
                idx = np.abs(s11_results['freq_ghz'] - proj_freq).argmin(); s11_at_proj = s11_results['s11_db'][idx]
                self.s11_ax.axvline(proj_freq, color='lime', linestyle='-', alpha=0.5)
                self.s11_ax.plot(proj_freq, s11_at_proj, 'lime', marker='o', markersize=8)
                self.s11_ax.annotate(f'Projeto: {proj_freq:.2f} GHz\nS11: {s11_at_proj:.1f} dB', xy=(proj_freq, s11_at_proj),
                                     xytext=(proj_freq + 0.05 * (freq_max - freq_min), s11_at_proj - 5),
                                     fontsize=8, color='lime', arrowprops=dict(arrowstyle='->', color='lime', alpha=0.7))
            
            bw_mask = s11_results['s11_db'] < -10
            if np.any(bw_mask):
                freq_band = s11_results['freq_ghz'][bw_mask]
                if len(freq_band) > 0:
                    bw_min, bw_max = freq_band.min(), freq_band.max()
                    bw_ghz = bw_max - bw_min
                    self.s11_ax.axvspan(bw_min, bw_max, alpha=0.2, color='green')
                    bw_percent = 100 * bw_ghz / proj_freq if proj_freq > 0 else 0
                    bw_info = f"Largura de Banda (S11 < -10 dB): {bw_ghz*1000:.1f} MHz ({bw_percent:.1f}%)\n"
                    bw_info += f"Faixa: {bw_min:.3f} GHz a {bw_max:.3f} GHz"
                    self.bw_info_label.config(text=bw_info)
                else: self.bw_info_label.config(text="Nenhum ponto atinge S11 < -10 dB")
            else: self.bw_info_label.config(text="Nenhum ponto atinge S11 < -10 dB")
            
            self.s11_ax.set_xlabel('Frequência (GHz)'); self.s11_ax.set_ylabel('S11 (dB)')
            self.s11_ax.set_title('Perda de Retorno vs. Frequência'); self.s11_ax.grid(True, alpha=0.3)
            self.s11_ax.set_xlim(freq_min, freq_max); self.s11_ax.set_ylim(min(s11_results['s11_db'].min(), -30), 0)
            self.s11_canvas.draw()
        except ValueError: Messagebox.show_error("Erro de Entrada", "Por favor, verifique os valores de frequência.")
        except Exception as e: Messagebox.show_error("Erro Inesperado", f"Ocorreu um erro: {e}")

    def load_demo(self):
        self.clear_fields()
        self.entries["freq"].insert(0, "2.4")
        self.entries["z0"].insert(0, "50")
        self.entries["zl_real"].insert(0, "75")
        self.entries["zl_imag"].insert(0, "25")
        self.entries["h"].insert(0, "1.6")
        self.entries["er"].insert(0, "4.4")
        self.calculate()

# =============================================================================
# BLOCO PRINCIPAL PARA INICIAR A APLICAÇÃO
# =============================================================================

if __name__ == "__main__":
    root = tb.Window(themename="superhero")
    app = CalculatorApp(root)
    root.mainloop()