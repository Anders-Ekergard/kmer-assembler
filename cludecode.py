"""
MVP: SNP-baserat avelsstod kopplat till foderkvot (FCR)

Steg:
1. Ladda in genotyp- och fenotypdata
2. Koppla ihop via djur_ID
3. Jamfor FCR mellan genotypgrupper per SNP-markor
4. Rakna avels-poang och rangordna djur
5. Visa resultat

Krav: pandas, scipy
"""

import pandas as pd
from scipy import stats


# ---------- 1. Ladda in data ----------

def ladda_genotyper(path):
    """Forvantat format (CSV): djur_id, markor1, markor2, ...
    Varje markor-varde ar t.ex. 'AA', 'AB' eller 'BB'."""
    return pd.read_csv(path)


def ladda_fenotyper(path):
    """Forvantat format (CSV): djur_id, fcr"""
    return pd.read_csv(path)


# ---------- 2. Koppla ihop data ----------

def koppla_data(genotyper, fenotyper):
    return genotyper.merge(fenotyper, on="djur_id", how="inner")


# ---------- 3. Jamfor FCR per genotypgrupp, per markor ----------

def analysera_markorer(data, markor_kolumner):
    resultat = []
    for markor in markor_kolumner:
        grupper = data.groupby(markor)["fcr"]

        # Hoppa over markorer med bara en genotyp representerad (gar inte att jamfora)
        if grupper.ngroups < 2:
            continue

        medelvarden = grupper.mean()
        basta_genotyp = medelvarden.idxmin()  # lagst FCR = mest foder-effektiv
        samplen = [grupp.values for _, grupp in grupper]

        # Enkel ANOVA (funkar bade for 2 och fler grupper)
        f_varde, p_varde = stats.f_oneway(*samplen)

        resultat.append({
            "markor": markor,
            "basta_genotyp": basta_genotyp,
            "medel_fcr_basta": round(medelvarden[basta_genotyp], 3),
            "p_varde": round(p_varde, 4),
            "signifikant": p_varde < 0.05,
        })

    return pd.DataFrame(resultat)


# ---------- 4. Rakna avels-poang och rangordna ----------

def berakna_poang(data, markor_analys, markor_kolumner):
    """Ett djur far +1 poang for varje SIGNIFIKANT markor dar det bar den basta genotypen."""
    signifikanta = markor_analys[markor_analys["signifikant"]]

    poang = pd.Series(0, index=data.index)
    for _, rad in signifikanta.iterrows():
        markor = rad["markor"]
        basta = rad["basta_genotyp"]
        poang += (data[markor] == basta).astype(int)

    data = data.copy()
    data["avels_poang"] = poang
    return data.sort_values(["avels_poang", "fcr"], ascending=[False, True])


# ---------- 5. Visa resultat ----------

def visa_resultat(markor_analys, rangordnade_djur, topp_n=5):
    print("== Markor-analys ==")
    if markor_analys.empty:
        print("Inga jamforbara markorer hittades.")
    else:
        print(markor_analys.to_string(index=False))

    print(f"\n== Topp {topp_n} avelskandidater ==")
    kolumner = ["djur_id", "fcr", "avels_poang"]
    print(rangordnade_djur[kolumner].head(topp_n).to_string(index=False))


# ---------- Korning ----------

if __name__ == "__main__":
    # Byt ut mot dina riktiga filer:
    # genotyper = ladda_genotyper("genotyper.csv")
    # fenotyper = ladda_fenotyper("fenotyper.csv")

    # -- Syntetiskt exempel sa du kan testa direkt --
    genotyper = pd.DataFrame({
        "djur_id": [f"R{i}" for i in range(1, 11)],
        "markor_A": ["AA", "AB", "BB", "AA", "AB", "BB", "AA", "AB", "BB", "AA"],
        "markor_B": ["AB", "AB", "AA", "BB", "AA", "AB", "AB", "BB", "AA", "AB"],
    })
    fenotyper = pd.DataFrame({
        "djur_id": [f"R{i}" for i in range(1, 11)],
        "fcr": [1.4, 1.6, 2.1, 1.3, 1.7, 2.0, 1.5, 1.8, 2.2, 1.35],
    })

    markor_kolumner = ["markor_A", "markor_B"]

    data = koppla_data(genotyper, fenotyper)
    markor_analys = analysera_markorer(data, markor_kolumner)
    rangordnade = berakna_poang(data, markor_analys, markor_kolumner)

    visa_resultat(markor_analys, rangordnade)