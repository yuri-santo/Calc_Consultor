# Calculadora SAP - Operações e Fórmulas (Python)

Aplicativo em Python (Tkinter) voltado para **conferência de cálculos** no contexto de processos tipo SAP:
- cálculos rápidos de impostos (modelos auditáveis),
- cenários financeiros (PRICE/SAC, juros compostos, juros com aportes),
- conversões de unidades (KG/UDS),
- compensação de débitos x créditos.

> Objetivo: servir como **ferramenta de apoio** para validar números e simulações.
> Não substitui a determinação fiscal completa do SAP (que depende de dados mestres, tabelas e legislação).

---

## Aviso importante (escopo realista)

Impostos e regras de negócio “de verdade” variam por:
- produto/serviço (NCM, item de serviço),
- CFOP/CST, regime tributário, benefícios, base reduzida, ST/FCP,
- UF origem/destino, operação (venda, devolução, exportação, importação),
- atualizações legais.

Este app implementa **modelos didáticos e auditáveis**, com entradas explícitas.
Se quiser “nível SAP”, o caminho é:
1) adicionar tabelas/config (CSV/JSON) de regras, alíquotas e exceções;  
2) criar determinadores por cenário (UF, CFOP, CST, NCM etc.);  
3) manter atualização legal.

---

## Requisitos
- Python 3.10+ (recomendado)
- Tkinter (normalmente já vem com o Python no Windows)

---

## Como rodar

```bat
pip install --upgrade pip
python calculadora.py
