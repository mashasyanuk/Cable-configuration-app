import tkinter as tk
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib
matplotlib.use('TkAgg')

# Константы
d_k = 0.009
m_k_kg = 0.35
W_npa_kg = 480
g = 9.80665
rho = 1025
#c_n = 1.2
#c_t = 0.02

def get_cn(Re):
    if Re < 1000:
        return 0.5
    elif Re < 10000:
        return 1.0 + 0.2 * (Re - 1000) / 9000
    else:
        return 1.2

def get_ct(Re):
    return 0.02 * (10000 / Re)**0.2

def calculate_configuration(v_x, L_k, N, P_x, P_y, P_z, K):
    R_x =  K * v_x**2
    W_npa = - W_npa_kg * 10
    l_k = L_k / N
    w_k = 2.79296

    phi = np.zeros(N)
    phi_y = np.zeros(N)
    phi_z = np.zeros(N)
    x = np.zeros(N + 1)
    y = np.zeros(N + 1)
    z = np.zeros(N + 1)

    Tx0 = R_x
    Ty0 = W_npa + P_y
    Tz0 = P_z
    T0 = np.sqrt(Tx0**2 + Ty0**2 + Tz0**2)

    cos_phi = np.clip(Tx0 / T0, -1, 1)
    phi[0] = np.arccos(cos_phi)
    
    cos_phi_y = np.clip(Ty0 / T0, -1, 1)
    phi_y[0] = np.arccos(cos_phi_y)

    cos_phi_z = np.clip(Tz0 / T0, -1, 1)
    phi_z[0] = np.arccos(cos_phi_z)

    Tx = Tx0
    Ty = Ty0
    Tz = Tz0
    T_next=T0

    for i in range(N-1):
        print(i)
        v_ni = v_x * np.sin(phi[i])
        v_ti = v_x * np.cos(phi[i])

        # Число Рейнольдса
        Re_n = abs(v_ni) * d_k / 1e-6
        Re_t = abs(v_ti) * d_k / 1e-6
    
    # Защита от нуля
        if Re_n < 1: Re_n = 1
        if Re_t < 1: Re_t = 1

        R_ni = 0.5 * get_cn(Re_n) * l_k * d_k * rho * v_ni**2
        R_ti = 0.5 *  get_ct(Re_t)* np.pi * l_k * d_k * rho * v_ti**2

        denom1 = np.sqrt(np.cos(phi_z[i])**2 + np.cos(phi_y[i])**2)

        cos_x = np.sqrt(np.cos(phi_y[i])**2 + np.cos(phi_z[i])**2)
        cos_y = -(np.cos(phi[i]) * np.cos(phi_y[i])) / denom1
        cos_z = -(np.cos(phi[i]) * np.cos(phi_z[i])) / denom1

        Tx_next = Tx + R_ti * np.cos(phi[i]) + R_ni * cos_x
        Ty_next = Ty -  (w_k)*l_k + R_ni * cos_y + R_ti * np.cos(phi_y[i])
        Tz_next = Tz + R_ti * np.cos(phi_z[i]) + R_ni * cos_z

        T_next = np.sqrt(Tx_next**2 + Ty_next**2 + Tz_next**2)

        cos_phi_next = np.clip(Tx_next / T_next, -1, 1)
        phi[i+1] = np.arccos(cos_phi_next)
        cos_phi_next_y = np.clip(Ty_next / T_next, -1, 1)
        phi_y[i+1] = np.arccos(cos_phi_next_y)
        cos_phi_next_z = np.clip(Tz_next / T_next, -1, 1)
        phi_z[i+1] = np.arccos(cos_phi_next_z)

        Tx, Ty, Tz = Tx_next, Ty_next, Tz_next

    x[N] = 0
    y[N] = 0
    z[N] = 0
    for i in range(N-1, -1, -1):
        x[i] = x[i+1] + l_k * np.cos(phi[i])
        y[i] = y[i+1] + l_k * np.cos(phi_y[i])
        z[i] = z[i+1] + l_k * np.cos(phi_z[i])

    # Инверсия Y 
    for i in range(N-1, -1, -1):
        x[i] *= -1
    
    if v_x==0:
        x[0]=0
        y[0]=-L_k
        z[0]=0
        coord_label_x.config(text=f"X = {x[0]:.4f} м")
        coord_label_y.config(text=f"Y = {y[0]:.4f} м")
        coord_label_z.config(text=f"Z = {z[0]:.4f} м")
    else:

        coord_label_x.config(text=f"X = {x[0]:.4f} м")
        coord_label_y.config(text=f"Y = {y[0]:.4f} м")
        coord_label_z.config(text=f"Z = {z[0]:.4f} м")

    return x, y, z, phi, phi_z, phi_y

def plot_graph():
    try:
        v_x = float(entry_vx.get())
        L_k = float(entry_Lk.get())
        N = int(entry_N.get())
        P_x = float(entry_Px.get())
        P_y = float(entry_Py.get())
        P_z = float(entry_Pz.get())
        K = float(entry_K.get())

        if v_x < 0 or L_k <= 0 or N <= 0:
            return

        x, y, z, phi, phi_z, phi_y = calculate_configuration(v_x, L_k, N, P_x, P_y, P_z, K)

        for widget in frame_graph.winfo_children():
            widget.destroy()

        fig = Figure(figsize=(8, 10))
        ax = fig.add_subplot(211)
        ax.plot(x, y, 'b-', linewidth=2)
        ax.scatter(0, 0, color='red', s=50, label='Судно')
        ax.scatter(x[0], y[0], color='green', s=50, label='НПА')
        ax.set_xlabel('x, м')
        ax.set_ylabel('y, м')
        ax.set_title('Плоскость XY')
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.axis('equal')

        ax2 = fig.add_subplot(212)
        ax2.plot(z, y, 'r-', linewidth=2)
        ax2.scatter(0, 0, color='red', s=50, label='Судно')
        ax2.scatter(z[0], y[0], color='green', s=50, label='НПА')
        ax2.set_xlabel('z, м')
        ax2.set_ylabel('y, м')
        ax2.set_title('Плоскость YZ')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        ax2.axis('equal')
        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=frame_graph)
        canvas.draw()
        canvas.get_tk_widget().pack()        

    except:
        pass

root = tk.Tk()
root.title("Расчет кабеля НПА")
root.geometry("1200x800")

frame_left = tk.Frame(root, padx=10, pady=10)
frame_left.pack(side=tk.LEFT, fill=tk.Y)

tk.Label(frame_left, text="Коэффициент K, м/с:").pack()
entry_K = tk.Entry(frame_left)
entry_K.insert(0, "0.5")
entry_K.pack()

tk.Label(frame_left, text="Скорость v_x, м/с:").pack()
entry_vx = tk.Entry(frame_left)
entry_vx.insert(0, "0.5")
entry_vx.pack()

tk.Label(frame_left, text="Длина L_k, м:").pack()
entry_Lk = tk.Entry(frame_left)
entry_Lk.insert(0, "1000")
entry_Lk.pack()

tk.Label(frame_left, text="Звеньев N:").pack()
entry_N = tk.Entry(frame_left)
entry_N.insert(0, "1000")
entry_N.pack()

tk.Label(frame_left, text="Сила P_x, Н:").pack()
entry_Px = tk.Entry(frame_left)
entry_Px.insert(0, "0")
entry_Px.pack()

tk.Label(frame_left, text="Сила P_y, Н:").pack()
entry_Py = tk.Entry(frame_left)
entry_Py.insert(0, "0")
entry_Py.pack()

tk.Label(frame_left, text="Сила P_z, Н:").pack()
entry_Pz = tk.Entry(frame_left)
entry_Pz.insert(0, "0")
entry_Pz.pack()


tk.Button(frame_left, text="Построить", command=plot_graph, bg="lightblue").pack(pady=20)

frame_coords = tk.Frame(frame_left, relief=tk.GROOVE, bd=2)
frame_coords.pack(pady=(0, 20), fill=tk.X, padx=5)
tk.Label(frame_coords, text="Координаты НПА", font=("Arial", 10, "bold")).pack(pady=(5, 10))
coord_label_x = tk.Label(frame_coords, text="X = --- м", font=("Arial", 9))
coord_label_x.pack()
coord_label_y = tk.Label(frame_coords, text="Y = --- м", font=("Arial", 9))
coord_label_y.pack()
coord_label_z = tk.Label(frame_coords, text="Z = --- м", font=("Arial", 9))
coord_label_z.pack()


frame_graph = tk.Frame(root, bg="white")
frame_graph.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)


root.mainloop()
