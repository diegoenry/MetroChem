# MetroChem
Reference code for MetroChem

O Produto mínimo viável do MetroChem será desenvolvido em 4 etapas principais ao longo de 2019->2020.

1 - MetroConf - Análise conformacional.  
2 - MetroDock - Docking Molecular.  
3 - MetroMD - Dinâmica Molecular.  
4 - MetroMC - Monte Carlo.  

## MetroConf 

Pré-proposta #1: Desenvolvimento de um Gerador de Conformações de mínimo de energia

Pré-proposta #2: Desenvolvimento de um Gerador de Conformações Bioativas
para compostos de interesse farmacológico.

### Motivação:
A geração de conformeros moleculares de baixa energia é uma tarefa chave
em diversas áreas da química computacional, modelagem molecular, e
quimioinformática. Eles são especialmente importantes no contexto da
química medicinal, onde a performance dos protocolos de
descobrimento/desenvolvimento de novos fármacos (docking,
shape-matching, pharmacophore-matching) depende largamente da capacidade
de geração de conformeros bioativos. Essa é a primeira das ferramentas,
em seguida iremos desenvolver uma de Docking Molecular que pode 
aproveitar desta 1a ferramenta.

### Motivação extra:
Apesar da análise conformacional ser um tema normal das aulas de química, o Brasil não possui o know-how de geração de conformeros para <b>moléculas bioativas</b>, e nenhuma ferramenta de software para tanto.

### Justificativa:
A busca por conformações de biomoléculas é uma área interdisciplinar.
Ela envolve o conhecimento biológico, químico, físico e necessita
pesadamente da das ciências de informação e computação. O
desenvolvimento de ferramentas nesta área envolve planejamento e
desenvolvimento de algoritmos de busca e de otimização, e abre espaço
para o desenvolvimento/aprimoramento de funções de pontuação baseadas em
física ou empíricas, incluindo aprendizado de máquina. O projeto também
abre o leque para novas abordagem para modelar os problemas, e
implementações em arquiteturas de computação distribuída, GPU.


Portanto, temos a oportunidade de desenvolver uma nova plataforma de
pesquisa de longo prazo em química computacional.


Eu sugiro uma colaboração Diego/Heder/Lobosco para o negócio ir para a
frente bem rápido com publicações com essa mentalidade recrutamento de
alunos para debater com alguma frequência as partes mais desafiadoras do
código e desenvolver as especificações e aplicações.

Eu desejo que os produtos sejam de código aberto e de livre uso
acadêmico, no entanto licenciavel para a indústria como produto ao
serviço. Se formos realmente bons, podemos gerar receita para o
departamento.


# Recomendações de artigo para o MetroConf.

## Estado da Arte.
Hawkins, P. C. D. (2017). Conformation Generation: The State of the Art.
Journal of Chemical Information and Modeling, 57(8), 1747–1756.
https://doi.org/10.1021/acs.jcim.7b00221

## Aplicação de High-Throughput
Lyu, J., Wang, S., Balius, T. E., Singh, I., Levit, A., Moroz, Y. S., ...
Irwin, J. J. (2019). Ultra-large library docking for discovering new
chemotypes. Nature. https://doi.org/10.1038/s41586-019-0917-9

## Conformações bioativas.
Gürsoy, O., & Smieško, M. (2017). Searching for bioactive conformations
of drug-like ligands with current force fields: how good are we? Journal
of Cheminformatics, 9(1), 29. https://doi.org/10.1186/s13321-017-0216-0

## OpenMM
http://openmm.org/
