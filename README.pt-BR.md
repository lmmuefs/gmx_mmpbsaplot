[🇺🇸 English](README.md) · **🇧🇷 Português**

# gmx_mmpbsa_plot

Gráficos de MM/PBSA para quem usa o **gmx_MMPBSA** (campos de força Amber).
Desenvolvido no **LMM (Laboratório de Modelagem Molecular), UEFS**.

O script varre as subpastas, detecta os arquivos de cada sistema e abre um menu para
escolher o que plotar. No fim grava um relatório `.txt` com as escolhas e os resultados.

## Requisitos

Python 3.8+ com `numpy`, `pandas` e `matplotlib`.

```bash
pip install numpy pandas matplotlib
```

## Estrutura esperada

Uma subpasta por sistema (réplica, ligante, etc.), com os arquivos de saída do `gmx_MMPBSA`:

```
meu_projeto/
├── sistema_A/
│   ├── FINAL_RESULTS_MMPBSA.csv   # energias por frame (gmx_MMPBSA ... -eo FINAL_RESULTS_MMPBSA.csv)
│   └── FINAL_DECOMP_MMPBSA.csv    # decomposição por resíduo (-deo, com &decomp no input)
│                                  # (ou FINAL_DECOMP_MMPBSA.dat, já com média e desvio)
├── sistema_B/
│   └── ...
└── gmx_mmpbsa_plot.py
```

As médias ± desvio são calculadas a partir das energias por frame (desvio entre frames). Nem toda
pasta precisa ter os dois arquivos: o menu só oferece o que existe. As energias já estão em
kcal/mol (unidade do gmx_MMPBSA). A coluna `Frame #` é o frame da trajetória; por isso o menu
pergunta **quantos ns vale um frame** (`--ns-per-frame`, ex.: `0.01` para 10 ps) para montar o eixo
de tempo. Se o CSV tiver GB e PB, escolha com `--model GB|PB`. O ligante (`L:` no gmx_MMPBSA) é
detectado automaticamente.

## Uso

```bash
python gmx_mmpbsa_plot.py                 # menu interativo, pasta atual
python gmx_mmpbsa_plot.py /caminho/dados  # outra pasta raiz
```

O primeiro menu pergunta o **idioma** (Português do Brasil ou English US). A escolha vale
para tudo: menus e mensagens, textos dos gráficos, cabeçalhos das tabelas e relatório
`.txt`. Em português, os números dos eixos usam vírgula decimal. O último idioma usado
fica salvo no `plot_style.ini`. Para pular a pergunta: `--lang pt` ou `--lang en`.
(O texto do `--help` é sempre em inglês.)

Depois o menu pergunta, nesta ordem:

1. **Sistemas** a incluir (`1,3`, `1-3` ou `all`).
2. **Legenda** de cada sistema. O padrão é o nome da pasta; Enter mantém.
3. **Modelo**: MM/GBSA (Generalized Born) ou MM/PBSA (Poisson-Boltzmann). O script lê qual
   está no CSV e já sugere; escolher um que não está no arquivo dá erro. Com GBSA, os gráficos
   incluem também **ΔG gás (GGAS)** e **ΔG solvatação (GSOLV)**.
4. **Gráficos**: energia no tempo (total, vdW, eletrostática, polar, apolar; no GBSA também
   gás e solvatação), médias ± desvio, decomposição por resíduo.
5. **Tempo de cada frame** em ns e **faixa de tempo** (veja "Sobre o eixo do tempo" abaixo).
6. **Resíduos**: Top 5, Top 10 ou Todos. *Todos não é recomendado*: com centenas de
   resíduos o gráfico de barras fica ilegível, então só se gera o perfil por posição
   na sequência e o CSV completo.
7. **Ligante**: o ligante é detectado pelo prefixo `L:` do gmx_MMPBSA; o menu pergunta se
   mantém no ranking. Mantido, ele costuma dominar o topo.
8. **Cores**: hex separados por vírgula (`#0072B2,#D55E00`) ou Enter para a paleta padrão.
9. **Organização** das pastas de saída (ver abaixo).
10. **Formatos** (png, pdf, svg), exportar CSV dos dados, e janela interativa.

### Sobre o eixo do tempo

O CSV do gmx_MMPBSA só tem o **número do frame da trajetória** (`Frame #`), sem tempo; por
isso o script pergunta quantos ns vale um frame da trajetória (`--ns-per-frame`).
Exemplo: `dt = 0.002` ps e `nstxout-compressed = 5000` dão 10 ps por frame, ou seja,
`0.01` ns; o frame 50000 vira 500 ns. É o tempo de um frame da trajetória *original*, não o
intervalo entre duas linhas do CSV (o `interval` já está nos números dos frames). A pergunta
só aparece quando algum gráfico de energia no tempo é escolhido.

### Saída

Por padrão, tudo vai para uma pasta nova `MM_PBSA_plots_AAAAMMDD_HHMMSS/`, criada a cada
execução (então nada é sobrescrito), organizada por tipo de gráfico:

```
MM_PBSA_plots_AAAAMMDD_HHMMSS/
├── energy_vs_time/     # energia em função do tempo (um gráfico por componente)
├── summary/            # médias ± desvio dos componentes
├── per_residue/        # decomposição por resíduo (um por sistema)
├── tables/             # CSVs: rankings, tabelas "média ± desvio", dados dos gráficos
└── GB_mmpbsa_report_AAAAMMDD_HHMMSS.txt
```

Os **nomes dos arquivos dentro** ganham o prefixo do método: `GB_` ou `PB_`
(ex.: `GB_energy_total_vs_time.png`, `GB_mmpbsa_report_....txt`). O nome da pasta não muda.

Use `-o nome_da_pasta` para escolher a pasta de saída (se ela já existir, é reaproveitada,
e arquivos com o mesmo nome são sobrescritos).

Em português as pastas ficam `energia_vs_tempo/`, `medias/`, `por_residuo/` e `tabelas/`.
Outras organizações disponíveis: **por formato** (`png/`, `pdf/`, `svg/`, útil na hora
de submeter) ou **tudo na mesma pasta**.

### O que é o gráfico de médias (`energy_components_mean`)

Barras com a média de cada termo da energia de ligação, por sistema, lidas do
`FINAL_RESULTS_MMPBSA.csv`: vdW, eletrostática, solvatação polar, solvatação apolar e total
(total = soma dos quatro). No **MM/GBSA** as barras incluem também **ΔG gás** (GGAS = vdW +
eletrostática) e **ΔG solvatação** (GSOLV = polar + apolar), como se costuma reportar em artigos de GBSA. A média é sobre **toda a trajetória** e a barra de erro é o
**desvio padrão entre frames** (a flutuação da energia), não a incerteza da média.
Isso precisa constar na legenda da figura. O termo entrópico (−TΔS) não está incluído.

## Editando legendas, rótulos e fonte

Na primeira execução é criado o `plot_style.ini` na pasta raiz. Edite e rode de novo:

```ini
[names]
# aceita mathtext do matplotlib (comentários só em linha própria)
sistema_A = $\it{Mtb}$PS - KS1

[text_pt]
energy_ylabel = ENERGIA (kcal/mol)
residue_title = Decomposição de energia por resíduo
comp_vdw = vDW

[font]
family = Arial
size = 12
```

Os textos ficam em dois blocos, `[text_en]` e `[text_pt]`, um por idioma.
Um `plot_style.ini` de versões anteriores é atualizado sozinho na primeira execução
(as legendas que você definiu são mantidas).
A fonte padrão é **Arial**. Se ela não estiver instalada (comum em Linux), o script usa
o Liberation Sans, que tem as mesmas métricas, e registra isso no relatório.

Para ajustes finos fora do script:

- **SVG e PDF** mantêm o texto editável. Abra no Inkscape ou Illustrator
  (no SVG a fonte precisa estar instalada na máquina de quem abre).
- **Janela interativa** (opção do menu): o botão *Figure options* da barra do
  matplotlib altera título, rótulos, limites e cores antes de salvar.
- **CSV** dos dados: para refazer o gráfico em Origin, Excel ou Grace.

## Uso sem menu

Qualquer opção pula a pergunta correspondente. Com `-y` nada é perguntado: ele usa os
padrões (todos os sistemas, todos os formatos, ligante removido), então `--plots` é
obrigatório, e `--time`, `--ns-per-frame` e `--top` também quando esses gráficos forem escolhidos.
`--model GB|PB` pula a pergunta do modelo.

```bash
python gmx_mmpbsa_plot.py . --systems 1,3 --plots total,vdw,residues \
    --time 50 500 --ns-per-frame 0.01 --model GB --top 10 --ligand drop --formats png,svg --csv \
    --lang pt --layout type -y
```

`python gmx_mmpbsa_plot.py --help` lista todas as opções.

## Observações

- **Variáveis por modelo.** Polar = `EGB` (GB) ou `EPB` (PB); apolar = `ESURF` (GB) ou
  `ENPOLAR` + `EDISPER` / `ECAVITY` + `EDISPER` (PB). O script avisa se o total diferir da soma
  dos quatro componentes (termos de ligação/1-4 não nulos no delta).
- **Decomposição**: além do `FINAL_DECOMP_MMPBSA.csv` (por frame), o script aceita o
  `FINAL_DECOMP_MMPBSA.dat` (média e desvio já agregados pelo gmx_MMPBSA). Se os dois existirem,
  usa o CSV. O cabeçalho do `.dat` informa se foi `idecomp = 2` (termos 1-4 somados ao vdW e à
  eletrostática); vale citar no artigo. Só aparecem os resíduos do `print_res`.

- O gráfico de médias usa o `FINAL_RESULTS_MMPBSA.csv`, ou seja, a **trajetória inteira**.
  As médias do relatório para a faixa de tempo escolhida são calculadas a partir do CSV por frame
  e podem diferir.
- **Decomposição por resíduo e faixa de tempo.** A faixa de tempo escolhida vale também para a
  decomposição por resíduo: com o `FINAL_DECOMP_MMPBSA.csv` (por frame), média e desvio de cada
  resíduo usam só os frames dentro da faixa (tempo = `Frame #` × tempo por frame, como nos
  gráficos de energia no tempo). O `FINAL_DECOMP_MMPBSA.dat` já vem agregado e não pode ser
  cortado por tempo: nesse caso usa-se a trajetória inteira e aparece um aviso. O relatório
  informa a faixa usada nos resíduos.
- A tabela `*_top_N_table.csv` traz os valores no formato `média ± desvio`.

## Autoria

- Késsia Souza Santos: LMM (Laboratório de Modelagem Molecular), UEFS
