# Calculadora de Casamento (Stub, LC) com Carta de Smith Interativa

Aplicação educacional e prática para casamento de impedâncias com:
- Carta de Smith interativa
- Cálculo de stub paralelo em curto
- Redes de casamento L (LC)
- Dimensionamento de linha de microfita
- Varredura de S11 (largura de banda)
- Visualização do layout e corte da microfita

## Requisitos
- Python 3.8+
- Dependências:
  - tkinter (já incluído no Python; no Linux pode exigir pacote do sistema)
  - ttkbootstrap
  - numpy
  - matplotlib

No Linux (Ubuntu/Debian), instale Tk caso necessário:
- sudo apt-get update && sudo apt-get install -y python3-tk

## Instalação
1) Crie e ative um ambiente virtual (recomendado):
- python -m venv .venv
- Windows: .venv\Scripts\activate
- Linux/macOS: source .venv/bin/activate
2) Instale as dependências:
- pip install ttkbootstrap numpy matplotlib

## Execução
- python magrf.py

## Como Usar

A interface possui 4 abas principais.

1) Casamento com Stub
- Entradas:
  - Frequência (GHz)
  - Z0 (Ω)
  - ZL (Real) (Ω)
  - ZL (Imag) (Ω)
  - (Opcional) Altura do substrato h (mm) e εr para microfita
- Botões:
  - Calcular Stub: executa o cálculo e plota na Carta de Smith
  - Limpar: limpa campos e gráficos
  - Demo: preenche valores de exemplo
- Saídas:
  - Distância d e comprimento l do stub (em cm e frações de λ)
  - Carta de Smith com:
    - Ponto da carga (ZL)
    - Círculo de VSWR
    - Rotação ao longo da linha (d)
    - Conversão para admitância
    - Efeito do stub até o centro (casamento)
  - Dica: passe o mouse sobre a carta para ver Γ, Z e Y no rodapé

2) Casamento LC
- Entradas: Frequência (GHz), Z0 (Ω), ZL (Real/Imag) (Ω)
- Ações:
  - Calcular Redes LC: lista todas as soluções L-match possíveis
  - Para cada solução, use “Plotar & Ver Solução” para visualizar:
    - Trajetória na Carta de Smith
    - Diagrama esquemático (série/shunt, C/L) com valores
- Observação: a ordem da topologia é “Componente próximo à carga, Componente próximo à fonte”

3) Análise de Largura de Banda
- Requer uma solução de stub calculada
- Entradas:
  - Frequência inicial/final (GHz) e passo (GHz)
- Saídas:
  - Curva S11(dB) vs. frequência
  - Linha de -10 dB e destaque da faixa que satisfaz S11 < -10 dB
  - Marca da frequência de projeto

4) Visualização da Microfita
- Mostra:
  - Corte transversal (substrato, faixa condutora) com h e W
  - Visão superior (linha principal, stub e curto)
- Requer:
  - h (mm) e εr preenchidos
  - Uma solução de stub válida (usa d e l)

## Unidades e Convenções
- Frequência: GHz
- Impedâncias: Ω
- d, l: exibidos em cm e frações de λ
- Microfita:
  - Entrada: h (mm), εr (adimensional)
  - Saída: W (mm)
- Carta de Smith usa impedâncias normalizadas por Z0

## Exemplo Rápido (Demo)
1) Clique em “Demo”
2) Clique em “Calcular Stub”
3) Explore a Carta de Smith e, se desejar, vá à aba “Análise de Largura de Banda” e gere o S11
4) Preencha h=1.6 e εr=4.4 para ver a visualização da microfita

## Notas de Modelagem
- Stub e linha: modelo ideal, sem perdas
- Microfita: aproximação clássica (ex.: Hammerstad-Jensen) para W; resultados são aproximações
- LC: elementos ideais (L, C), sem parasitas; podem existir múltiplas soluções para o mesmo caso

## Solução de Problemas
- Erro de backend Tk/TkAgg ou janela não abre:
  - Instale o Tk (Linux): sudo apt-get install python3-tk
- Gráfico não atualiza ou erro de conversão:
  - Verifique se todos os campos numéricos estão preenchidos (use ponto “.” como separador decimal)
- S11 não plota:
  - Calcule o stub primeiro na aba principal

## Créditos
- tkinter + ttkbootstrap (UI)
- numpy (cálculos)
- matplotlib (Carta de Smith e gráficos)

Sugestões e melhorias são bem-vindas.

---

## Compilação/Empacotamento (Linux)
Este projeto inclui um Makefile para facilitar:

- Preparar ambiente e instalar dependências: `make install`
- Executar a aplicação: `make run`
- Gerar binário standalone (PyInstaller): `make build`
  - Saída: `dist/calculadora-casamento` (modo janela/GUI)

Requisitos extras para build: `python3-venv`, `python3-tk`, e o PyInstaller (instalado automaticamente pelo alvo `build`).

## Publicação no GitHub (passos sugeridos)
1. Iniciar repositório e primeiro commit:
   - git init
   - git add .
   - git commit -m "chore: initial commit"
   - git branch -M main
2. Crie um repositório no GitHub (web) e copie a URL (ex.: `https://github.com/<user>/<repo>.git`)
3. Vincule e envie:
   - git remote add origin <URL>
   - git push -u origin main
