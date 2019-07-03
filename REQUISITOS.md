# Especificação de requisitos
O presente documento tem o objetivo de especificar e estabelecer os requisitos para o desenvolvimento 
do programa MetroConf que tem como finalidade a análise conformacional de moléculas.
O primeiro método a ser implementado é a a geração de conformeros pela varredura dos ângulos dihedrais, a partir de uma geometria inicial.


## Descrição molecular;
Pela mecânica clássica, a descrição de uma molécula pode ser aproximada por um sistema de massas e molas a qual atribuimos parâmetros.
As massas são os átomos que possuem como parâmetros:
1. Tipo
2. Carga

Cada <i>tipo</i> de átomo definirá todos os demais parâmetros atômicos, e para as <b>molas</b> abaixo.
1. Epsilon
2. Sigma

As molas representadas pelas ligações químicas e interações não ligadas.  
<b>Interações ligadas</b>
1. Ligações
2. Ângulos
3. Dihedros próprios (torções)
4. Dihedros impróprios (fora do plano)

<b>Interaçes não ligadas</b>  
1. Forças de van der Waals.
2. Forças de Coulomb.

## A - Dados de entrada  
1. Leitura do arquivo de topologia  (butane.psf)  
2. Leitura de Coordenadas Arquivo.PDB  (butane.pdb)  
3. Leitura do arquivo de parâmetros (butane.frcmod)  


B- Output
B1 - Escrever a trajetoria da simulação.
B2 - Escrever a energia de um passo da simulação
B3 - Escrever o sumario das energias (média) ao final da simulação


C - Listas de interações ligadas (lista estática)
C1 - Gerar a lista de ligações ao combinar a topologia com os parâmetros de ligação (Bonds)
C2 - Gerar a lista de ligações ao combinar a topologia com os parâmetros de ângulo (Angles)
C3 - Gerar a lista de ligações ao combinar a topologia com os parâmetros de ângulos dihedros proprios (Torções)
C4 - Gerar a lista de ligações ao combinar a topologia com os parâmetros de ângulo dihedros impróprios (impropers, "out of plane")
C5 - Gerar a lista de vizinhos excluidos para cada átomo para NÃO calcular as interações de longo alcance. Os vizinhos excluidos são todos envolvidos em interações ligadas.
C6 - Gerar o potencial de vdW para cada TIPO de par de átomos usando regra combinatória.

D - Listas de interações não ligadas (lista dinâmica, atualizada a cada passo da simulação)
Na primeira implementação, fazer com distance based, e o quanto antes evoluir para uma busca GRID-Based.
D1 - Calcular a distância entre todos pares de átomos.
D2 - Gerar lista de vizinhos de acordo com um raio de corte para as interações de van der Waals..
D3 - Gerar lista de vizinhos de acordo com um raio de corte para as interações eletrostáticas.
*D4 - Para sistemas solvatados, aplicar condições periódicas de contorno e determinar o vizinho mais próximo pelo método da imagem mínima.


E - Interações ligadas (bonds, angles, dihedrals, impropers).
E1 - Calcular a energia e forças em cada átomo, devido às ligações (bonds)
E2 - Calcular a energia e forças em cada átomo, devido aos ângulos.
E3 - Calcular a energia e forças em cada átomo, devido aos dihedros proprios (torções)
E4 - Calcular a energia e forças em cada átomo, devido aos dihedros próprios (impropers)

F -  Interações de longo alcance ( Não ligadas)
F1 - Calcular a energia e forças em cada átomo, devido às interações de van der Waals.
F2 - Calcular a energia e forças em cada átomo, devido às interações eletrostáticas, pelo potencial de coulomb.
F3 - Calcular a energia e forças em cada átomo, devido às interações eletrostáticas, pelo potencial de coulomb com campo de reação.
F4 - Calcular a energia e forças em cada átomo, devido às interações eletrostáticas, pelo Particle Mesh Ewald.

G - Otimizar o sistema, com um minimizador a energia
G1 - Algoritmo simples de Steepest Descent.
G2 - Gradientes conjugados.
G3 - L-BFGS  https://en.wikipedia.org/wiki/Limited-memory_BFGS

H - Buscar conformação
H1.1 - Sistemático
H1.2 - Sistemático com minimização
H2.1 - AG
H2.2 - AG com formação nicho
