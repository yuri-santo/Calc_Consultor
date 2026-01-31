import tkinter as tk
from tkinter import ttk
import ast
import math


# =========================
# Helpers de entrada/saída
# =========================
def to_float(txt: str) -> float:
    """
    Aceita:
      - 1234,56
      - 1.234,56
      - 1234.56
    """
    s = (txt or "").strip()
    if not s:
        raise ValueError("campo vazio")

    # caso BR com milhar + vírgula decimal
    if "," in s:
        s = s.replace(".", "")
        s = s.replace(",", ".")
    return float(s)


def to_int(txt: str) -> int:
    v = to_float(txt)
    return int(round(v))


def money_br(x: float) -> str:
    s = f"{x:,.2f}"  # 1,234.56
    s = s.replace(",", "X").replace(".", ",").replace("X", ".")
    return f"R$ {s}"


def pct(x: float) -> str:
    # x em percentual, ex: 12.5 -> "12,5%"
    s = f"{x:.6f}".rstrip("0").rstrip(".")
    s = s.replace(".", ",")
    return f"{s}%"


def num(x: float) -> str:
    if abs(x - round(x)) < 1e-12:
        return str(int(round(x)))
    return f"{x:.12f}".rstrip("0").rstrip(".")


def r2(x: float) -> float:
    return round(float(x), 2)


# =========================
# Expressão segura (sem eval solto)
# =========================
FUNCOES = {
    "sqrt": math.sqrt,
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "log": math.log,     # log(x) ou log(x, base)
    "log10": math.log10,
    "abs": abs,
    "round": round,
}

CONSTANTES = {"pi": math.pi, "e": math.e}
OPS_BIN = (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.Mod)
OPS_UN = (ast.UAdd, ast.USub)


def _deg_to_rad_if_needed(func: str, args: list[float], usar_graus: bool) -> list[float]:
    if usar_graus and func in ("sin", "cos", "tan") and args:
        args = args[:]
        args[0] = math.radians(args[0])
    return args


def avaliar_expressao(expr: str, usar_graus: bool = True) -> float:
    expr = (expr or "").strip()
    if not expr:
        raise ValueError("expressão vazia")

    expr = expr.replace(",", ".")
    tree = ast.parse(expr, mode="eval")

    def walk(n):
        if isinstance(n, ast.Expression):
            return walk(n.body)

        if isinstance(n, ast.Constant):
            if isinstance(n.value, (int, float)):
                return float(n.value)
            raise ValueError("constante inválida")

        if isinstance(n, ast.Name):
            if n.id in CONSTANTES:
                return float(CONSTANTES[n.id])
            raise ValueError(f"nome não permitido: {n.id}")

        if isinstance(n, ast.UnaryOp):
            if not isinstance(n.op, OPS_UN):
                raise ValueError("sinal não permitido")
            v = walk(n.operand)
            return +v if isinstance(n.op, ast.UAdd) else -v

        if isinstance(n, ast.BinOp):
            if not isinstance(n.op, OPS_BIN):
                raise ValueError("operador não permitido")

            a = walk(n.left)
            b = walk(n.right)

            if isinstance(n.op, ast.Add):
                return a + b
            if isinstance(n.op, ast.Sub):
                return a - b
            if isinstance(n.op, ast.Mult):
                return a * b
            if isinstance(n.op, ast.Div):
                return a / b
            if isinstance(n.op, ast.Pow):
                return a ** b
            if isinstance(n.op, ast.Mod):
                return a % b

            raise ValueError("operação inválida")

        if isinstance(n, ast.Call):
            if not isinstance(n.func, ast.Name):
                raise ValueError("chamada inválida")

            nome = n.func.id
            if nome not in FUNCOES:
                raise ValueError(f"função não permitida: {nome}")

            args = [walk(a) for a in n.args]
            args = _deg_to_rad_if_needed(nome, args, usar_graus)

            return float(FUNCOES[nome](*args))

        raise ValueError("expressão não suportada")

    return walk(tree)


# =========================
# Motor base de imposto (genérico)
# =========================
def imposto_base_padrao(total_b: float, exclusao: float, outros: float, red_base_pct: float) -> float:
    """
    BaseAmt = (TotalB - Exclusão - Outros) * (1 - RedBase)
    com trava de não ficar negativo.
    """
    if red_base_pct < 0 or red_base_pct > 100:
        raise ValueError("redução de base fora do intervalo (0 a 100)")
    base = (total_b - exclusao - outros) * (1 - red_base_pct / 100.0)
    return r2(max(0.0, base))


def imposto_valor(base_amt: float, rate_pct: float) -> float:
    if rate_pct < 0:
        raise ValueError("alíquota inválida")
    return r2(base_amt * (rate_pct / 100.0))


def imposto_final(tax_amt: float, credit: float) -> float:
    return r2(tax_amt - credit)


# =========================
# Operações (cada uma retorna lista de (label, valor_formatado))
# =========================

# --- Básicos ---
def op_percentual(p):
    v = p["Valor"]
    r = p["Percentual (%)"]
    res = v * (r / 100.0)
    return [("Valor", money_br(v)), ("Percentual", pct(r)), ("Resultado", money_br(r2(res)))]


def op_somar_taxa(p):
    base = p["Valor base"]
    rate = p["Taxa/Alíquota (%)"]
    tax = base * (rate / 100.0)
    total = base + tax
    return [
        ("Base", money_br(base)),
        ("Taxa", pct(rate)),
        ("Taxa (valor)", money_br(r2(tax))),
        ("Total", money_br(r2(total))),
    ]


def op_desconto(p):
    base = p["Valor"]
    rate = p["Desconto (%)"]
    d = base * (rate / 100.0)
    total = base - d
    return [
        ("Valor", money_br(base)),
        ("Desconto", pct(rate)),
        ("Desconto (valor)", money_br(r2(d))),
        ("Total", money_br(r2(total))),
    ]


# --- Conversões / logística ---
def op_kg_para_uds(p):
    kg = p["Quantidade (kg)"]
    kg_por_ud = p["Kg por unidade"]
    if kg_por_ud <= 0:
        raise ValueError("kg por unidade precisa ser > 0")
    uds = kg / kg_por_ud
    return [
        ("Kg", num(kg)),
        ("Kg por unidade", num(kg_por_ud)),
        ("Unidades (UDS)", num(uds)),
    ]


def op_uds_para_kg(p):
    uds = p["Quantidade (UDS)"]
    kg_por_ud = p["Kg por unidade"]
    if kg_por_ud <= 0:
        raise ValueError("kg por unidade precisa ser > 0")
    kg = uds * kg_por_ud
    return [
        ("UDS", num(uds)),
        ("Kg por unidade", num(kg_por_ud)),
        ("Kg total", num(kg)),
    ]


def op_compensacao_credito(p):
    debito = p["Débitos (R$)"]
    credito = p["Créditos (R$)"]
    saldo = r2(debito - credito)
    a_pagar = r2(max(0.0, saldo))
    a_compensar = r2(max(0.0, -saldo))
    return [
        ("Débitos", money_br(debito)),
        ("Créditos", money_br(credito)),
        ("Saldo (débito - crédito)", money_br(saldo)),
        ("A pagar", money_br(a_pagar)),
        ("Crédito a compensar", money_br(a_compensar)),
    ]


# --- Financeiro ---
def op_juros_simples(p):
    pv = p["Principal (PV)"]
    i = p["Taxa ao mês (%)"]
    n = to_int(str(p["Meses"]))
    juros = pv * (i / 100.0) * n
    fv = pv + juros
    return [
        ("PV", money_br(pv)),
        ("Taxa", pct(i) + " a.m."),
        ("Meses", str(n)),
        ("Juros", money_br(r2(juros))),
        ("FV", money_br(r2(fv))),
    ]


def op_juros_compostos(p):
    pv = p["Principal (PV)"]
    i = p["Taxa ao mês (%)"] / 100.0
    n = to_int(str(p["Meses"]))
    fv = pv * ((1 + i) ** n)
    juros = fv - pv
    return [
        ("PV", money_br(pv)),
        ("Taxa", pct(i * 100) + " a.m."),
        ("Meses", str(n)),
        ("Juros", money_br(r2(juros))),
        ("FV", money_br(r2(fv))),
    ]


def op_juros_compostos_aportes(p):
    """
    FV com aporte mensal (PMT). Por padrão, aporte no fim do mês.

    FV = PV*(1+i)^n + PMT * [((1+i)^n - 1)/i]
    """
    pv = p["Valor inicial (PV)"]
    pmt = p["Aporte mensal (PMT)"]
    i = p["Taxa ao mês (%)"] / 100.0
    n = to_int(str(p["Meses"]))

    if n < 0:
        raise ValueError("meses inválido")
    if abs(i) < 1e-12:
        fv = pv + pmt * n
    else:
        fv = pv * ((1 + i) ** n) + pmt * (((1 + i) ** n - 1) / i)

    total_aportado = pv + pmt * n
    ganho = fv - total_aportado

    return [
        ("PV", money_br(pv)),
        ("PMT (aporte)", money_br(pmt)),
        ("Taxa", pct(i * 100) + " a.m."),
        ("Meses", str(n)),
        ("Total aportado", money_br(r2(total_aportado))),
        ("FV", money_br(r2(fv))),
        ("Ganho (FV - aportes)", money_br(r2(ganho))),
    ]


def op_price(p):
    pv = p["Valor financiado (PV)"]
    i = p["Taxa ao mês (%)"] / 100.0
    n = to_int(str(p["Nº parcelas"]))

    if n <= 0:
        raise ValueError("nº parcelas inválido")

    if abs(i) < 1e-12:
        pmt = pv / n
    else:
        pmt = pv * (i * (1 + i) ** n) / ((1 + i) ** n - 1)

    total_pago = pmt * n
    juros_total = total_pago - pv

    return [
        ("PV", money_br(pv)),
        ("Taxa", pct(i * 100) + " a.m."),
        ("Parcelas", str(n)),
        ("PMT (parcela)", money_br(r2(pmt))),
        ("Total pago", money_br(r2(total_pago))),
        ("Juros total", money_br(r2(juros_total))),
    ]


def op_sac(p):
    pv = p["Valor financiado (PV)"]
    i = p["Taxa ao mês (%)"] / 100.0
    n = to_int(str(p["Nº parcelas"]))

    if n <= 0:
        raise ValueError("nº parcelas inválido")

    amort = pv / n
    saldo = pv
    juros_total = 0.0
    primeira = ultima = 0.0

    for k in range(1, n + 1):
        juros = saldo * i
        prest = amort + juros
        juros_total += juros

        if k == 1:
            primeira = prest
        if k == n:
            ultima = prest

        saldo -= amort

    total_pago = pv + juros_total

    return [
        ("PV", money_br(pv)),
        ("Taxa", pct(i * 100) + " a.m."),
        ("Parcelas", str(n)),
        ("Amortização", money_br(r2(amort))),
        ("1ª parcela", money_br(r2(primeira))),
        ("Última parcela", money_br(r2(ultima))),
        ("Juros total", money_br(r2(juros_total))),
        ("Total pago", money_br(r2(total_pago))),
    ]


# --- Brasil (modelo simplificado e auditável) ---
def op_brasil_icms_proprio(p):
    # BaseICMS = ValorOperacao - Exclusoes
    # BaseFinal = BaseICMS * (1 - RedBase%)
    # ICMS = BaseFinal * Aliquota
    valor = p["Valor da operação (R$)"]
    exclusoes = p["Exclusões (R$)"]
    red = p["Redução de base (%)"]
    aliq = p["Alíquota ICMS (%)"]

    base = imposto_base_padrao(valor, exclusoes, 0.0, red)
    icms = imposto_valor(base, aliq)
    return [
        ("Valor operação", money_br(valor)),
        ("Exclusões", money_br(exclusoes)),
        ("Redução base", pct(red)),
        ("Base ICMS", money_br(base)),
        ("Alíquota", pct(aliq)),
        ("ICMS", money_br(icms)),
    ]


def op_brasil_icms_difal(p):
    # DIFAL (simplificado):
    # ICMS_dest = Base * Aliq_dest
    # ICMS_orig = Base * Aliq_orig
    # DIFAL = ICMS_dest - ICMS_orig
    base = p["Base ICMS (R$)"]
    aliq_dest = p["Alíquota destino (%)"]
    aliq_orig = p["Alíquota origem (%)"]

    icms_dest = imposto_valor(base, aliq_dest)
    icms_orig = imposto_valor(base, aliq_orig)
    difal = r2(icms_dest - icms_orig)

    return [
        ("Base", money_br(base)),
        ("ICMS destino", money_br(icms_dest)),
        ("ICMS origem", money_br(icms_orig)),
        ("DIFAL", money_br(difal)),
    ]


def op_brasil_icms_st(p):
    """
    Modelo didático:
      BaseST = (ValorOperacao * (1 + MVA%)) + IPI - Desconto
      ICMS_ST = BaseST * Aliq_ST
      ICMS_Proprio = (ValorOperacao - Exclusoes) * Aliq_Propria  (aqui simplificado)
      ICMS_A_Recolher = ICMS_ST - ICMS_Proprio
    """
    valor = p["Valor da operação (R$)"]
    mva = p["MVA/Markup (%)"]
    ipi = p["IPI (R$)"]
    desconto = p["Desconto (R$)"]
    aliq_st = p["Alíquota ST (%)"]
    aliq_prop = p["Alíquota própria (%)"]
    exclusoes = p["Exclusões (R$)"]

    base_st = r2((valor * (1 + mva / 100.0)) + ipi - desconto)
    icms_st = imposto_valor(base_st, aliq_st)

    base_prop = r2(max(0.0, valor - exclusoes))
    icms_prop = imposto_valor(base_prop, aliq_prop)

    recolher = r2(icms_st - icms_prop)

    return [
        ("Valor operação", money_br(valor)),
        ("MVA", pct(mva)),
        ("IPI", money_br(ipi)),
        ("Desconto", money_br(desconto)),
        ("Base ST", money_br(base_st)),
        ("ICMS ST", money_br(icms_st)),
        ("Base própria", money_br(base_prop)),
        ("ICMS próprio", money_br(icms_prop)),
        ("ICMS a recolher (ST - próprio)", money_br(recolher)),
    ]


def op_brasil_ipi(p):
    base = p["Base IPI (R$)"]
    aliq = p["Alíquota IPI (%)"]
    red = p["Redução base (%)"]
    base_f = imposto_base_padrao(base, 0.0, 0.0, red)
    ipi = imposto_valor(base_f, aliq)
    total = r2(base + ipi)
    return [
        ("Base original", money_br(base)),
        ("Redução", pct(red)),
        ("Base final", money_br(base_f)),
        ("IPI", money_br(ipi)),
        ("Total (base + IPI)", money_br(total)),
    ]


def op_brasil_iss(p):
    valor = p["Valor do serviço (R$)"]
    exclusoes = p["Exclusões legais (R$)"]
    aliq = p["Alíquota ISS (%)"]
    base = r2(max(0.0, valor - exclusoes))
    iss = imposto_valor(base, aliq)
    return [
        ("Valor serviço", money_br(valor)),
        ("Exclusões", money_br(exclusoes)),
        ("Base ISS", money_br(base)),
        ("ISS", money_br(iss)),
    ]


def op_brasil_pis(p):
    """
    PIS:
      Cumulativo: PIS = Base * 1,65%
      Não-cumulativo: PIS_deb = Base * 7,65%; PIS_dev = deb - credito
    """
    regime = int(p["Regime (1=cumulativo, 2=nao cumulativo)"])
    base = p["Base (R$)"]
    credito = p["Créditos (R$)"]

    if regime == 1:
        aliq = 1.65
        pis = imposto_valor(base, aliq)
        return [
            ("Regime", "Cumulativo"),
            ("Base", money_br(base)),
            ("Alíquota", pct(aliq)),
            ("PIS", money_br(pis)),
        ]

    if regime == 2:
        aliq = 7.65
        deb = imposto_valor(base, aliq)
        devido = r2(deb - credito)
        return [
            ("Regime", "Não-cumulativo"),
            ("Base", money_br(base)),
            ("Débito PIS", money_br(deb)),
            ("Créditos", money_br(credito)),
            ("PIS devido", money_br(devido)),
        ]

    raise ValueError("regime inválido (use 1 ou 2)")


def op_brasil_cofins(p):
    """
    COFINS:
      Cumulativa: 3,0% (varia na prática; aqui deixei campo editável)
      Não-cumulativa: 7,6% (com créditos)
    """
    regime = int(p["Regime (1=cumulativo, 2=nao cumulativo)"])
    base = p["Base (R$)"]
    credito = p["Créditos (R$)"]

    if regime == 1:
        aliq = p["Alíquota cumulativa (%)"]
        cof = imposto_valor(base, aliq)
        return [
            ("Regime", "Cumulativo"),
            ("Base", money_br(base)),
            ("Alíquota", pct(aliq)),
            ("COFINS", money_br(cof)),
        ]

    if regime == 2:
        aliq = 7.60
        deb = imposto_valor(base, aliq)
        devido = r2(deb - credito)
        return [
            ("Regime", "Não-cumulativo"),
            ("Base", money_br(base)),
            ("Débito COFINS", money_br(deb)),
            ("Créditos", money_br(credito)),
            ("COFINS devida", money_br(devido)),
        ]

    raise ValueError("regime inválido (use 1 ou 2)")


def op_brasil_irpj(p):
    """
    IRPJ (simplificado):
      Lucro real: IRPJ = lucro * 15% + adicional 10% sobre excedente (20.000)
      Lucro presumido: Base = receita * presuncao%; IRPJ = base*15% + adicional
    """
    modo = int(p["Modo (1=lucro real, 2=presumido)"])

    if modo == 1:
        lucro = p["Lucro operacional (R$)"]
        ir = r2(lucro * 0.15)
        adicional = r2(max(0.0, lucro - 20000.0) * 0.10)
        total = r2(ir + adicional)
        return [
            ("Modo", "Lucro real"),
            ("Lucro", money_br(lucro)),
            ("IRPJ (15%)", money_br(ir)),
            ("Adicional (10% excedente)", money_br(adicional)),
            ("IRPJ total", money_br(total)),
        ]

    if modo == 2:
        receita = p["Receita bruta (R$)"]
        pres = p["% presunção (%)"] / 100.0
        base = r2(receita * pres)
        ir = r2(base * 0.15)
        adicional = r2(max(0.0, base - 20000.0) * 0.10)
        total = r2(ir + adicional)
        return [
            ("Modo", "Lucro presumido"),
            ("Receita", money_br(receita)),
            ("Presunção", pct(p["% presunção (%)"])),
            ("Base presumida", money_br(base)),
            ("IRPJ (15%)", money_br(ir)),
            ("Adicional", money_br(adicional)),
            ("IRPJ total", money_br(total)),
        ]

    raise ValueError("modo inválido (use 1 ou 2)")


def op_brasil_csll(p):
    lucro = p["Base CSLL/Lucro (R$)"]
    aliq = p["Alíquota CSLL (%)"]
    csll = imposto_valor(lucro, aliq)
    return [
        ("Base", money_br(lucro)),
        ("Alíquota", pct(aliq)),
        ("CSLL", money_br(csll)),
    ]


# --- Global / importação (modelo simplificado) ---
def op_vat_gst(p):
    base = p["Base (R$)"]
    rate = p["Taxa VAT/GST (%)"]
    tax = imposto_valor(base, rate)
    total = r2(base + tax)
    return [("Base", money_br(base)), ("Taxa", pct(rate)), ("VAT/GST", money_br(tax)), ("Total", money_br(total))]


def op_sales_tax(p):
    base = p["Base (R$)"]
    rate = p["Taxa Sales Tax (%)"]
    tax = imposto_valor(base, rate)
    total = r2(base + tax)
    return [("Base", money_br(base)), ("Taxa", pct(rate)), ("Sales Tax", money_br(tax)), ("Total", money_br(total))]


def op_customs_import(p):
    """
    Importação (didático):
      CIF = mercadoria + frete + seguro
      II = CIF * tarifa_ii
      IPI = (CIF + II) * aliq_ipi
      ICMS = (CIF + II + IPI) * aliq_icms   (simplificado, sem "por dentro")
      PIS = (CIF + II + IPI) * aliq_pis
      COFINS = (CIF + II + IPI) * aliq_cofins
    """
    merc = p["Valor mercadoria (R$)"]
    frete = p["Frete (R$)"]
    seguro = p["Seguro (R$)"]
    ii_rate = p["Tarifa II (%)"]
    ipi_rate = p["Alíquota IPI (%)"]
    icms_rate = p["Alíquota ICMS (%)"]
    pis_rate = p["Taxa PIS (%)"]
    cof_rate = p["Taxa COFINS (%)"]

    cif = r2(merc + frete + seguro)
    ii = imposto_valor(cif, ii_rate)
    ipi = imposto_valor(r2(cif + ii), ipi_rate)
    base_trib = r2(cif + ii + ipi)

    icms = imposto_valor(base_trib, icms_rate)
    pis = imposto_valor(base_trib, pis_rate)
    cof = imposto_valor(base_trib, cof_rate)

    total_impostos = r2(ii + ipi + icms + pis + cof)
    total = r2(cif + total_impostos)

    return [
        ("CIF", money_br(cif)),
        ("II", money_br(ii)),
        ("IPI", money_br(ipi)),
        ("Base (CIF+II+IPI)", money_br(base_trib)),
        ("ICMS", money_br(icms)),
        ("PIS", money_br(pis)),
        ("COFINS", money_br(cof)),
        ("Total impostos", money_br(total_impostos)),
        ("Total importação", money_br(total)),
    ]


# --- Cascata “SAP style” (simplificada, mas com passos) ---
def op_cascata_nf(p):
    """
    Sequência didática parecida com o que você descreveu:
      Base ICMS -> ICMS
      IPI sobre base
      PIS/COFINS sobre (base + IPI) (simplificado)
      ISS opcional
      Total = base + impostos
    """
    base_nf = p["Valor base NF (R$)"]
    excl = p["Exclusões (R$)"]
    desc = p["Descontos incond. (R$)"]
    red_icms = p["Redução base ICMS (%)"]
    aliq_icms = p["Alíquota ICMS (%)"]
    aliq_ipi = p["Alíquota IPI (%)"]
    regime = int(p["Regime PIS/COFINS (1=cum, 2=nao cum)"])
    credit_pis = p["Crédito PIS (R$)"]
    credit_cof = p["Crédito COFINS (R$)"]
    aliq_iss = p["Alíquota ISS (%)"]
    valor_serv = p["Valor serviço (R$)"]

    # Base ICMS
    base_icms_bruta = r2(max(0.0, base_nf - excl - desc))
    base_icms = imposto_base_padrao(base_icms_bruta, 0.0, 0.0, red_icms)
    icms = imposto_valor(base_icms, aliq_icms)

    # IPI (sobre base_nf)
    ipi = imposto_valor(base_nf, aliq_ipi)

    # PIS/COFINS sobre (base_nf + ipi) simplificado
    base_pc = r2(base_nf + ipi)
    if regime == 1:
        pis = imposto_valor(base_pc, 1.65)
        cof = imposto_valor(base_pc, 7.60)
    elif regime == 2:
        pis_deb = imposto_valor(base_pc, 7.65)
        cof_deb = imposto_valor(base_pc, 7.60)
        pis = r2(pis_deb - credit_pis)
        cof = r2(cof_deb - credit_cof)
    else:
        raise ValueError("regime inválido (1 ou 2)")

    # ISS opcional (se valor_serv > 0)
    iss_base = r2(max(0.0, valor_serv))
    iss = imposto_valor(iss_base, aliq_iss) if iss_base > 0 else 0.0

    total_impostos = r2(icms + ipi + pis + cof + iss)
    total = r2(base_nf + total_impostos)

    out = [
        ("Base NF", money_br(base_nf)),
        ("Base ICMS bruta", money_br(base_icms_bruta)),
        ("Base ICMS final", money_br(base_icms)),
        ("ICMS", money_br(icms)),
        ("IPI", money_br(ipi)),
        ("Base PIS/COFINS (base+IPI)", money_br(base_pc)),
        ("PIS", money_br(pis)),
        ("COFINS", money_br(cof)),
        ("ISS", money_br(iss)),
        ("Total impostos", money_br(total_impostos)),
        ("Total final", money_br(total)),
    ]

    return out


# =========================
# Catálogo de operações (UI)
# =========================
OPERACOES = [
    # Expressões (fica numa aba separada)
    # Básicos
    {"cat": "Básico", "nome": "Calcular X% de um valor", "campos": [("Valor", "money"), ("Percentual (%)", "float")], "fn": op_percentual},
    {"cat": "Básico", "nome": "Somar taxa/alíquota (base + %)", "campos": [("Valor base", "money"), ("Taxa/Alíquota (%)", "float")], "fn": op_somar_taxa},
    {"cat": "Básico", "nome": "Aplicar desconto (valor - %)", "campos": [("Valor", "money"), ("Desconto (%)", "float")], "fn": op_desconto},

    # Conversões e logística
    {"cat": "Logística / Unidades", "nome": "Converter KG -> UDS", "campos": [("Quantidade (kg)", "float"), ("Kg por unidade", "float")], "fn": op_kg_para_uds},
    {"cat": "Logística / Unidades", "nome": "Converter UDS -> KG", "campos": [("Quantidade (UDS)", "float"), ("Kg por unidade", "float")], "fn": op_uds_para_kg},

    # Compensação
    {"cat": "Apuração / Compensação", "nome": "Compensar débitos x créditos", "campos": [("Débitos (R$)", "money"), ("Créditos (R$)", "money")], "fn": op_compensacao_credito},

    # Financeiro
    {"cat": "Financeiro", "nome": "Juros simples", "campos": [("Principal (PV)", "money"), ("Taxa ao mês (%)", "float"), ("Meses", "int")], "fn": op_juros_simples},
    {"cat": "Financeiro", "nome": "Juros compostos", "campos": [("Principal (PV)", "money"), ("Taxa ao mês (%)", "float"), ("Meses", "int")], "fn": op_juros_compostos},
    {"cat": "Financeiro", "nome": "Juros compostos + aportes mensais", "campos": [("Valor inicial (PV)", "money"), ("Aporte mensal (PMT)", "money"), ("Taxa ao mês (%)", "float"), ("Meses", "int")], "fn": op_juros_compostos_aportes},
    {"cat": "Financeiro", "nome": "Financiamento PRICE (parcela fixa)", "campos": [("Valor financiado (PV)", "money"), ("Taxa ao mês (%)", "float"), ("Nº parcelas", "int")], "fn": op_price},
    {"cat": "Financeiro", "nome": "Financiamento SAC (parcela decrescente)", "campos": [("Valor financiado (PV)", "money"), ("Taxa ao mês (%)", "float"), ("Nº parcelas", "int")], "fn": op_sac},

    # Brasil
    {"cat": "Brasil - ICMS/IPI/ISS", "nome": "ICMS próprio (com exclusões e redução)", "campos": [("Valor da operação (R$)", "money"), ("Exclusões (R$)", "money"), ("Redução de base (%)", "float"), ("Alíquota ICMS (%)", "float")], "fn": op_brasil_icms_proprio},
    {"cat": "Brasil - ICMS/IPI/ISS", "nome": "ICMS DIFAL (destino - origem)", "campos": [("Base ICMS (R$)", "money"), ("Alíquota destino (%)", "float"), ("Alíquota origem (%)", "float")], "fn": op_brasil_icms_difal},
    {"cat": "Brasil - ICMS/IPI/ISS", "nome": "ICMS-ST (com MVA e IPI)", "campos": [
        ("Valor da operação (R$)", "money"),
        ("MVA/Markup (%)", "float"),
        ("IPI (R$)", "money"),
        ("Desconto (R$)", "money"),
        ("Alíquota ST (%)", "float"),
        ("Alíquota própria (%)", "float"),
        ("Exclusões (R$)", "money"),
    ], "fn": op_brasil_icms_st},
    {"cat": "Brasil - ICMS/IPI/ISS", "nome": "IPI (normal com redução opcional)", "campos": [("Base IPI (R$)", "money"), ("Alíquota IPI (%)", "float"), ("Redução base (%)", "float")], "fn": op_brasil_ipi},
    {"cat": "Brasil - ICMS/IPI/ISS", "nome": "ISS (serviços)", "campos": [("Valor do serviço (R$)", "money"), ("Exclusões legais (R$)", "money"), ("Alíquota ISS (%)", "float")], "fn": op_brasil_iss},

    {"cat": "Brasil - PIS/COFINS", "nome": "PIS (cumulativo ou não-cumulativo)", "campos": [("Regime (1=cumulativo, 2=nao cumulativo)", "int"), ("Base (R$)", "money"), ("Créditos (R$)", "money")], "fn": op_brasil_pis},
    {"cat": "Brasil - PIS/COFINS", "nome": "COFINS (cumulativo ou não-cumulativo)", "campos": [
        ("Regime (1=cumulativo, 2=nao cumulativo)", "int"),
        ("Base (R$)", "money"),
        ("Créditos (R$)", "money"),
        ("Alíquota cumulativa (%)", "float"),
    ], "fn": op_brasil_cofins},

    {"cat": "Brasil - Lucro", "nome": "IRPJ (lucro real ou presumido)", "campos": [
        ("Modo (1=lucro real, 2=presumido)", "int"),
        ("Lucro operacional (R$)", "money"),
        ("Receita bruta (R$)", "money"),
        ("% presunção (%)", "float"),
    ], "fn": op_brasil_irpj},
    {"cat": "Brasil - Lucro", "nome": "CSLL (alíquota configurável)", "campos": [("Base CSLL/Lucro (R$)", "money"), ("Alíquota CSLL (%)", "float")], "fn": op_brasil_csll},

    # Global / Import
    {"cat": "Global", "nome": "VAT/GST padrão", "campos": [("Base (R$)", "money"), ("Taxa VAT/GST (%)", "float")], "fn": op_vat_gst},
    {"cat": "Global", "nome": "Sales Tax (EUA - simples)", "campos": [("Base (R$)", "money"), ("Taxa Sales Tax (%)", "float")], "fn": op_sales_tax},
    {"cat": "Importação", "nome": "Importação (CIF + II + IPI + ICMS + PIS + COFINS)", "campos": [
        ("Valor mercadoria (R$)", "money"),
        ("Frete (R$)", "money"),
        ("Seguro (R$)", "money"),
        ("Tarifa II (%)", "float"),
        ("Alíquota IPI (%)", "float"),
        ("Alíquota ICMS (%)", "float"),
        ("Taxa PIS (%)", "float"),
        ("Taxa COFINS (%)", "float"),
    ], "fn": op_customs_import},

    # Cascata
    {"cat": "SAP / Cascata", "nome": "Cascata NF (ICMS + IPI + PIS/COFINS + ISS)", "campos": [
        ("Valor base NF (R$)", "money"),
        ("Exclusões (R$)", "money"),
        ("Descontos incond. (R$)", "money"),
        ("Redução base ICMS (%)", "float"),
        ("Alíquota ICMS (%)", "float"),
        ("Alíquota IPI (%)", "float"),
        ("Regime PIS/COFINS (1=cum, 2=nao cum)", "int"),
        ("Crédito PIS (R$)", "money"),
        ("Crédito COFINS (R$)", "money"),
        ("Alíquota ISS (%)", "float"),
        ("Valor serviço (R$)", "money"),
    ], "fn": op_cascata_nf},
]


# =========================
# UI
# =========================
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Calculadora SAP - Operações e Fórmulas")
        self.geometry("640x560")
        self.minsize(620, 540)

        self.usar_graus = tk.BooleanVar(value=True)
        self.inputs = {}

        self._build()

    def _build(self):
        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=10, pady=10)

        tab_expr = ttk.Frame(nb)
        tab_ops = ttk.Frame(nb)

        nb.add(tab_expr, text="Calculadora")
        nb.add(tab_ops, text="Operações SAP")

        self._build_expr(tab_expr)
        self._build_ops(tab_ops)

    # --- Aba calculadora por expressão ---
    def _build_expr(self, parent):
        top = ttk.Frame(parent)
        top.pack(fill="x", pady=(5, 10))

        self.inp_expr = ttk.Entry(top, font=("Segoe UI", 14))
        self.inp_expr.pack(fill="x")
        self.inp_expr.focus_set()

        opts = ttk.Frame(parent)
        opts.pack(fill="x", pady=(0, 8))
        ttk.Checkbutton(opts, text="Trigonometria em graus", variable=self.usar_graus).pack(side="left")

        self.lbl_expr = ttk.Label(parent, text="Resultado: ", font=("Segoe UI", 11))
        self.lbl_expr.pack(fill="x", pady=(0, 10))

        grid = ttk.Frame(parent)
        grid.pack(fill="both", expand=True)

        for c in range(4):
            grid.columnconfigure(c, weight=1)

        buttons = [
            ("7", "7"), ("8", "8"), ("9", "9"), ("/", "/"),
            ("4", "4"), ("5", "5"), ("6", "6"), ("*", "*"),
            ("1", "1"), ("2", "2"), ("3", "3"), ("-", "-"),
            ("0", "0"), (".", "."), ("+", "+"), ("**", "**"),
            ("(", "("), (")", ")"), ("pi", "pi"), ("e", "e"),
            ("sqrt(", "sqrt("), ("sin(", "sin("), ("cos(", "cos("), ("tan(", "tan("),
            ("log(", "log("), ("abs(", "abs("), ("round(", "round("), ("%", "%"),
        ]

        r = c = 0
        for txt, ins in buttons:
            ttk.Button(grid, text=txt, command=lambda s=ins: self._ins_expr(s)).grid(
                row=r, column=c, sticky="nsew", padx=4, pady=4, ipady=6
            )
            c += 1
            if c == 4:
                c = 0
                r += 1

        actions = ttk.Frame(parent)
        actions.pack(fill="x", pady=(10, 0))

        ttk.Button(actions, text="Calcular (Enter)", command=self.calc_expr).pack(side="left", expand=True, fill="x", padx=(0, 5))
        ttk.Button(actions, text="Limpar (Esc)", command=self.clear_expr).pack(side="left", expand=True, fill="x", padx=(5, 0))

        self.bind("<Return>", lambda _: self.calc_expr())
        self.bind("<Escape>", lambda _: self.clear_expr())

    def _ins_expr(self, s: str):
        self.inp_expr.insert(tk.END, s)
        self.inp_expr.focus_set()

    def clear_expr(self):
        self.inp_expr.delete(0, tk.END)
        self.lbl_expr.config(text="Resultado: ")
        self.inp_expr.focus_set()

    def calc_expr(self):
        try:
            val = avaliar_expressao(self.inp_expr.get(), usar_graus=self.usar_graus.get())
            self.lbl_expr.config(text=f"Resultado: {num(val)}")
        except Exception as e:
            self.lbl_expr.config(text=f"Resultado: erro ({e})")

    # --- Aba operações ---
    def _build_ops(self, parent):
        header = ttk.Frame(parent)
        header.pack(fill="x", pady=(5, 8))

        ttk.Label(header, text="Categoria:").pack(side="left")
        self.cmb_cat = ttk.Combobox(header, state="readonly", width=26)
        self.cmb_cat.pack(side="left", padx=(6, 14))

        ttk.Label(header, text="Operação:").pack(side="left")
        self.cmb_op = ttk.Combobox(header, state="readonly")
        self.cmb_op.pack(side="left", fill="x", expand=True, padx=(6, 0))

        cats = sorted({o["cat"] for o in OPERACOES})
        self.cmb_cat["values"] = cats
        self.cmb_cat.current(0)
        self.cmb_cat.bind("<<ComboboxSelected>>", lambda _: self._load_ops_for_cat())

        self.cmb_op.bind("<<ComboboxSelected>>", lambda _: self._render_fields())

        # corpo em duas colunas: inputs + resultado
        body = ttk.Frame(parent)
        body.pack(fill="both", expand=True)

        left = ttk.Frame(body)
        left.pack(side="left", fill="both", expand=False, padx=(0, 8))

        right = ttk.Frame(body)
        right.pack(side="left", fill="both", expand=True)

        # campos dinâmicos
        self.frm_fields = ttk.LabelFrame(left, text="Entradas")
        self.frm_fields.pack(fill="both", expand=True)

        btns = ttk.Frame(left)
        btns.pack(fill="x", pady=(8, 0))
        ttk.Button(btns, text="Calcular", command=self.calc_op).pack(side="left", expand=True, fill="x", padx=(0, 4))
        ttk.Button(btns, text="Limpar", command=self.clear_op).pack(side="left", expand=True, fill="x", padx=(4, 0))

        # saída
        self.out = tk.Text(right, wrap="word")
        self.out.pack(fill="both", expand=True)

        self._log("Dica: aceita 1.234,56 e 1234,56. Campos inteiros aceitam 12 ou 12,0.\n")
        self._load_ops_for_cat()

    def _log(self, msg: str):
        self.out.insert(tk.END, msg)
        self.out.see(tk.END)

    def _load_ops_for_cat(self):
        cat = self.cmb_cat.get()
        ops = [o for o in OPERACOES if o["cat"] == cat]
        self.ops_current = ops
        self.cmb_op["values"] = [o["nome"] for o in ops]
        self.cmb_op.current(0)
        self._render_fields()

    def _render_fields(self):
        for w in self.frm_fields.winfo_children():
            w.destroy()
        self.inputs = {}

        op = self.ops_current[self.cmb_op.current()]
        campos = op["campos"]

        self.frm_fields.columnconfigure(1, weight=1)

        for i, (nome, tipo) in enumerate(campos):
            ttk.Label(self.frm_fields, text=nome + ":").grid(row=i, column=0, sticky="w", padx=8, pady=6)
            ent = ttk.Entry(self.frm_fields, width=28)
            ent.grid(row=i, column=1, sticky="ew", padx=8, pady=6)
            self.inputs[nome] = (ent, tipo)

        # foco no primeiro campo
        if campos:
            first = campos[0][0]
            self.inputs[first][0].focus_set()

    def clear_op(self):
        for ent, _ in self.inputs.values():
            ent.delete(0, tk.END)

    def calc_op(self):
        op = self.ops_current[self.cmb_op.current()]
        fn = op["fn"]

        try:
            params = {}
            for nome, (ent, tipo) in self.inputs.items():
                raw = ent.get()
                if tipo == "int":
                    params[nome] = to_int(raw)
                else:
                    params[nome] = to_float(raw)

            res = fn(params)

            self._log(f"\n[{op['cat']} -> {op['nome']}]\n")
            for k, v in res:
                self._log(f"{k}: {v}\n")

        except Exception as e:
            self._log(f"\nErro: {e}\n")


if __name__ == "__main__":
    App().mainloop()
