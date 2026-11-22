#!/usr/bin/env python3
"""
QUPLAS - bilancio di flusso semplificato: sorgente e+ -> Ps- -> fotodistacco -> Ps neutro
==========================================================================================

Parametri presi da: Galanti et al., Radiation Physics and Chemistry 251 (2027) 114388
("Design of an experimental setup for continuous positronium beam formation via Ps-
production and laser photo-detachment"), salvo quelli marcati  # IPOTESI.

Ogni valore marcato  # IPOTESI  NON e' nel paper: e' una stima di ordine di grandezza
da sostituire con i valori reali del gruppo. Il modello e' volutamente semplice
(1D/2D, non relativistico dove trascurabile) e serve a confrontare scenari,
non a sostituire le simulazioni del gruppo (Vinelli et al. 2023, Vinelli 2024).

Uso:   python3 quplas_bilancio_flusso.py
Modificare i dizionari SCENARI e PARAMETRI qui sotto per esplorare.
"""

import math

# ---------------------------------------------------------------------------
# Costanti fisiche
# ---------------------------------------------------------------------------
C = 2.99792458e8            # m/s
ME_C2 = 510998.95           # eV, energia a riposo dell'elettrone
HC_EV_NM = 1239.84198       # eV*nm
EV = 1.602176634e-19        # J

# ---------------------------------------------------------------------------
# Parametri dal paper
# ---------------------------------------------------------------------------
TAU_PSM = 479e-12           # s, vita media Ps- (Ceeh et al. 2011)
TAU_OPS = 142e-9            # s, vita media orto-Ps 1^3S1
TAU_2S = 1.14e-6            # s, vita di annichilazione 2^3S1
F_ORTO = 0.75               # frazione orto-Ps dopo fotodistacco (1:3)
ETA_CONV_MAX = 0.015        # efficienza max e+ -> Ps- (Terabe et al. 2012)
AFFINITA_PS = 0.32688       # eV, affinita' elettronica del Ps (Michishio et al. 2020)
SIGMA_PD = 6.5e-17 * 1e-4   # m^2, sezione d'urto fotodistacco vicino al picco
LAMBDA = 1560e-9            # m, laser
P_CIRC = 197.9e3            # W, potenza circolante (simulazione "cavita' calda")
L_CAV = 2.04679             # m, lunghezza cavita' misurata
R_MIRROR = 2.0              # m, raggio di curvatura nominale
TARGET_PS = 5e7             # Ps/s per 10% su g in 11 mesi (Vinelli et al. 2023)
MESI_RIF = 11

# ---------------------------------------------------------------------------
# Scenari di sorgente
# ---------------------------------------------------------------------------
ATT_50MCI = 50e-3 * 3.7e10            # Bq
BETA_PLUS_NA22 = 0.90                 # rapporto di ramificazione beta+ del 22Na

SCENARI = {
    # nome: positroni moderati al secondo
    "22Na 50 mCi + moderatore W":   ATT_50MCI * BETA_PLUS_NA22 * 2e-4,   # IPOTESI eff. W ~1e-4..5e-4
    "22Na 50 mCi + neon solido":    ATT_50MCI * BETA_PLUS_NA22 * 5e-3,   # ~0.5% (letteratura)
    "NEPOMUC (Monaco)":             9e8,                                 # fascio primario 1 keV
    "Sorgente futura 1e10":         1e10,                                # limite alto citato nel paper
}

# ---------------------------------------------------------------------------
# Parametri di progetto (modificabili)
# ---------------------------------------------------------------------------
PARAMETRI = dict(
    eta_trasporto=0.5,      # IPOTESI: trasporto + focalizzazione e+ sul convertitore
    eta_conv=ETA_CONV_MAX,  # e+ -> Ps- (1.5% = miglior valore di letteratura)
    E_psm_eV=3000.0,        # IPOTESI: energia cinetica Ps- dopo accelerazione
    gap_mm=1.0,             # IPOTESI: gap convertitore-griglia (campo uniforme)
    drift_mm=3.0,           # IPOTESI: griglia -> asse del modo laser (>= ~3 w di distanza)
    angolo_incrocio_deg=90, # angolo tra velocita' Ps- e asse cavita' (90 = perpendicolare)
    sigma_spot_mm=0.5,      # IPOTESI: rms del fascio e+ sul convertitore (per asse)
    E_trasv_emiss_eV=0.2,   # IPOTESI: energia trasversa rms di emissione Ps- (per asse)
    # --- catena estesa, tutti IPOTESI ---
    L_ecc_m=0.05,           # distanza fotodistacco -> eccitazione UV/MW
    eta_ecc=0.3,            # efficienza 1^3S -> 2^3S
    L_interf_m=0.5,         # lunghezza interferometro
    accettanza_mrad=2.0,    # semi-angolo di accettanza dell'interferometro
)

# ---------------------------------------------------------------------------
# Modello
# ---------------------------------------------------------------------------
def waist_cavita(L=L_CAV, R=R_MIRROR, lam=LAMBDA):
    """Waist del TEM00 in una cavita' simmetrica a due specchi curvi."""
    return math.sqrt(lam / math.pi * math.sqrt(L * (2 * R - L)) / 2)


def velocita(E_eV, massa_eV):
    gamma = 1 + E_eV / massa_eV
    beta = math.sqrt(1 - 1 / gamma**2)
    return beta * C, beta, gamma


def sopravvivenza_psm(E_eV, gap_mm, drift_mm):
    """Ps- parte (quasi) da fermo, accelera uniformemente nel gap, poi deriva fino al laser."""
    v, _, _ = velocita(E_eV, 3 * ME_C2)
    t_gap = 2 * gap_mm * 1e-3 / v
    t_drift = drift_mm * 1e-3 / v
    t = t_gap + t_drift
    return math.exp(-t / TAU_PSM), t


def prob_fotodistacco(E_eV, angolo_deg, sigma_fascio_m, w=None, P=P_CIRC, n_b=201, n_phi=64):
    """
    Probabilita' media di fotodistacco per un fascio Ps- gaussiano (rms sigma_fascio_m nella
    direzione perpendicolare al piano di incrocio) che attraversa un modo TEM00 in onda stazionaria.
    - perpendicolare: lo ione resta a fase fissa dell'onda stazionaria -> media su cos^2
    - obliquo: attraversa molte frange -> intensita' media = 2 x onda progressiva
    """
    w = w or waist_cavita()
    v, _, _ = velocita(E_eV, 3 * ME_C2)
    s = math.sin(math.radians(angolo_deg))
    h_nu = HC_EV_NM / (LAMBDA * 1e9) * EV
    # fluenza di fotoni (una direzione) per parametro d'urto b = 0
    I0 = 2 * P / (math.pi * w**2)
    F0 = I0 * w * math.sqrt(math.pi / 2) / (v * s) / h_nu
    perpendicolare = abs(math.cos(math.radians(angolo_deg))) * w * 3 < LAMBDA * 5

    tot, norm = 0.0, 0.0
    bmax = 5 * max(sigma_fascio_m, 1e-9)
    for i in range(n_b):
        b = -bmax + 2 * bmax * i / (n_b - 1)
        peso = math.exp(-b**2 / (2 * sigma_fascio_m**2)) if sigma_fascio_m > 0 else (1.0 if i == n_b // 2 else 0.0)
        F = F0 * math.exp(-2 * b**2 / w**2)
        if perpendicolare:
            p = sum(1 - math.exp(-SIGMA_PD * 4 * F * math.cos(math.pi * k / n_phi) ** 2)
                    for k in range(n_phi)) / n_phi
        else:
            p = 1 - math.exp(-SIGMA_PD * 2 * F)
        tot += peso * p
        norm += peso
    return tot / norm


def divergenze(E_eV, E_trasv_eV):
    """Angoli rms per asse: emissione dal convertitore e rinculo del fotodistacco (mrad)."""
    # emissione: v_trasv / v_long (angolo compresso dall'accelerazione)
    th_emiss = math.sqrt(E_trasv_eV / E_eV)
    # rinculo: Ps- -> Ps + e, energia in eccesso divisa 2:1 tra e e Ps
    E_foton = HC_EV_NM / (LAMBDA * 1e9)
    eccesso = E_foton - AFFINITA_PS
    p_rinculo = math.sqrt(4 / 3 * ME_C2 * eccesso)           # eV/c
    _, beta, gamma = velocita(E_eV, 3 * ME_C2)
    p_ps = 2 * ME_C2 * beta * gamma                           # eV/c, momento del Ps neutro
    th_rinc_max = p_rinculo / p_ps
    th_rinc_rms = th_rinc_max / math.sqrt(3)                  # distribuzione isotropa (semplificazione)
    return th_emiss * 1e3, th_rinc_rms * 1e3, th_rinc_max * 1e3


def bilancio(sorgente_eplus, p=PARAMETRI, verbose=False):
    E = p["E_psm_eV"]
    v, _, _ = velocita(E, 3 * ME_C2)
    n_conv = sorgente_eplus * p["eta_trasporto"]
    n_psm = n_conv * p["eta_conv"]
    surv, t = sopravvivenza_psm(E, p["gap_mm"], p["drift_mm"])
    n_psm_laser = n_psm * surv

    # dimensione del fascio Ps- all'incrocio
    v_t = math.sqrt(2 * p["E_trasv_emiss_eV"] / (3 * ME_C2)) * C
    sig = math.hypot(p["sigma_spot_mm"] * 1e-3, v_t * t)
    p_pd = prob_fotodistacco(E, p["angolo_incrocio_deg"], sig)
    n_ps = n_psm_laser * p_pd
    n_ops = n_ps * F_ORTO

    # catena estesa (IPOTESI)
    th_em, th_rec, _ = divergenze(E, p["E_trasv_emiss_eV"])
    th_tot = math.hypot(th_em, th_rec)
    f_acc = 1 - math.exp(-p["accettanza_mrad"] ** 2 / (2 * th_tot ** 2))
    v_ps = v  # la velocita' non cambia apprezzabilmente nel fotodistacco
    s_ops = math.exp(-p["L_ecc_m"] / (v_ps * TAU_OPS))
    s_2s = math.exp(-p["L_interf_m"] / (v_ps * TAU_2S))
    n_interf = n_ops * s_ops * p["eta_ecc"] * f_acc * s_2s

    r = dict(e_plus=sorgente_eplus, su_convertitore=n_conv, ps_meno=n_psm,
             ps_meno_al_laser=n_psm_laser, ps=n_ps, orto_ps=n_ops, in_interferometro=n_interf,
             sopravvivenza=surv, p_pd=p_pd, sigma_fascio_mm=sig * 1e3,
             th_emiss=th_em, th_rinc=th_rec, f_acc=f_acc)
    return r


def precisione_g(flusso, mesi=MESI_RIF):
    """Scala statistica 1/sqrt(N) rispetto al riferimento del paper (10% con 5e7 Ps/s in 11 mesi)."""
    return 0.10 * math.sqrt(TARGET_PS / flusso) * math.sqrt(MESI_RIF / mesi)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
def fmt(x):
    return f"{x:9.2e}"


def main():
    w = waist_cavita()
    print("=" * 78)
    print("QUPLAS - bilancio di flusso (modello semplificato)")
    print("=" * 78)
    print(f"Waist del modo in cavita': {w*1e3:.2f} mm  |  potenza circolante {P_CIRC/1e3:.0f} kW")
    p = PARAMETRI
    print(f"Ps- a {p['E_psm_eV']/1e3:.1f} keV, gap {p['gap_mm']} mm, drift {p['drift_mm']} mm, "
          f"incrocio {p['angolo_incrocio_deg']} gradi")

    print("\n1) CATENA PER SCENARIO (flussi al secondo)")
    print("-" * 78)
    hdr = f"{'scenario':30s}{'e+ mod.':>10s}{'Ps-':>10s}{'Ps- laser':>10s}{'o-Ps':>10s}{'x target':>10s}"
    print(hdr)
    for nome, s in SCENARI.items():
        r = bilancio(s)
        print(f"{nome:30s}{fmt(r['e_plus'])} {fmt(r['ps_meno'])} {fmt(r['ps_meno_al_laser'])} "
              f"{fmt(r['orto_ps'])} {r['orto_ps']/TARGET_PS:9.1e}")
    r = bilancio(1.0)
    print(f"\nFattori: trasporto {p['eta_trasporto']:.2f} x conversione {p['eta_conv']:.3f} x "
          f"sopravvivenza Ps- {r['sopravvivenza']:.2f} x fotodistacco {r['p_pd']:.2f} x orto {F_ORTO}")
    print(f"=> {r['orto_ps']:.2e} o-Ps per positrone moderato")
    print(f"Fascio Ps- all'incrocio: sigma = {r['sigma_fascio_mm']:.2f} mm (waist laser {w*1e3:.2f} mm)")

    print("\n2) PRECISIONE SU g (scala 1/sqrt(N), confronto con o-Ps dopo fotodistacco)")
    print("   ATTENZIONE: valida solo se il target 5e7 Ps/s del paper e' definito allo stesso")
    print("   punto della catena e alla stessa velocita' del fascio.")
    print("-" * 78)
    for nome, s in SCENARI.items():
        r = bilancio(s)
        print(f"{nome:30s} o-Ps {r['orto_ps']:9.2e}/s  ->  dg/g ~ {precisione_g(r['orto_ps'])*100:8.0f}% in 11 mesi")

    print("\n3) SCANSIONE IN ENERGIA DEI Ps- (per positrone moderato)")
    print("-" * 78)
    print(f"{'E [keV]':>8s}{'sopravv.':>10s}{'P_fotod.':>10s}{'prodotto':>10s}{'th_emiss':>10s}{'th_rinc':>10s}  [mrad rms]")
    for E in [0.5, 1, 2, 3, 5, 8, 12, 20]:
        q = dict(PARAMETRI, E_psm_eV=E * 1e3)
        r = bilancio(1.0, q)
        print(f"{E:8.1f}{r['sopravvivenza']:10.2f}{r['p_pd']:10.2f}{r['sopravvivenza']*r['p_pd']:10.2f}"
              f"{r['th_emiss']:10.1f}{r['th_rinc']:10.1f}")
    print("   (la sensibilita' dell'interferometro scala come T^2 ~ 1/v^2: energie alte la riducono)")

    print("\n4) ANGOLO DI INCROCIO FASCIO-LASER")
    print("-" * 78)
    for ang in [90, 45, 20, 10, 5]:
        q = dict(PARAMETRI, angolo_incrocio_deg=ang)
        r = bilancio(1.0, q)
        print(f"  {ang:3d} gradi: P_fotodistacco = {r['p_pd']:.2f}")

    print("\n5) DISTANZA GRIGLIA -> LASER")
    print("-" * 78)
    for d in [1, 2, 3, 5, 10, 20]:
        q = dict(PARAMETRI, drift_mm=d)
        r = bilancio(1.0, q)
        print(f"  {d:3d} mm: sopravvivenza Ps- = {r['sopravvivenza']:.2f}")

    print("\n6) CATENA ESTESA FINO ALL'INTERFEROMETRO (parametri quasi tutti IPOTESI)")
    print("-" * 78)
    r = bilancio(1.0)
    print(f"  accettanza {p['accettanza_mrad']} mrad su divergenza totale "
          f"{math.hypot(r['th_emiss'], r['th_rinc']):.1f} mrad rms -> frazione {r['f_acc']:.3f}")
    for nome, s in SCENARI.items():
        r = bilancio(s)
        print(f"  {nome:30s} {r['in_interferometro']:9.2e} Ps(2^3S)/s nell'interferometro")
    print("\nFine. Sostituire i valori # IPOTESI con quelli reali per un bilancio affidabile.")


if __name__ == "__main__":
    main()
