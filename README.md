README
========================================================

#### 2014-07-06

Passos dados:

 - Algoritmo para simulação dos fatores com volatilidade estocástica implementado;
 - Simulação das volatilidades estocásticas validada;
 - Simulação das cargas fatoriais na forma PLT para FSV programada corretamente e de maneira eficiente;
 - Amostragem dos parâmetros de autorregressão da volatilidade estocástica validada via Gibbs no arquivo `SampleSvTest.R`;
 - Testes feitos no arquivo `SampleFsvFactorsTest.R` para validar funções criadas.
 - No arquivo `tests/src/SampleSvItMn.cpp` foi testada a velocidade do algoritmo FFBS para simulação das volatilidades quando as variáveis latentes auxiliares da aproximação por mistura de normais é feita a cada iteração do filtro progressivo. Quando as variáveis latentes auxiliares são simuladas antes do filtro progressivo (usando a função `GetLogX2ApproxParms`), o algoritmo fica três vezes mais rápido do que quando se o aplica dentro do filtro.
 
Passos futuros:
 
 - Rever documentação das funções;
 - Construir FFBS para fatores dinâmicos com volatilidade estocástica;
 - Construir funções que calculam a verossimilhança marginal dos modelos FSV e DFSV para comparação de modelos;
 - Aplicar modelo de regressão fatorial;
 - Aplicar modelo de regressão fatorial com volatilidade estocástica (FSVR);
 - Aplicar modelo fatorial dinâmico;
 - Programar função que amostra cargas fatoriais a partir das restrições de Bai e Wang (2012) para aplicar modelo fatorial dinâmico multinível.
 - Aplicar modelo fatorial dinâmico multinível;
 - Aplicar modelo fatorial dinâmico de volatilidade estocástica;
 - Aplicar modelo fatorial dinâmico multinível com volatilidade estocástica;
 - Gerar dados artificiais para ilustração dos modelos;
 - Realizar análise de resíduos (correlogramas e covariâncias)
 
#### 2014-07-03

Passos dados:

 - Código da função de perda validado. A estimação dos parâmetros fica bastante próxima dos verdadeiros parâmetros. Há somente lentidão de convergência quando a ordem dinâmica dos fatores na equação das observações é próxima da ordem do VAR. Entretanto, as distribuições *a posteriori* são bastante próximas dos verdadeiros parâmetros.
 - Iniciada junção das funções que simulam a volatilidade estocástica e respectivos parâmetros para estimação do modelo fatorial de volatilidade estocástica.
 
Próximos passos:
  
  - Transcrever funções para tirar amostras da distribuição *a posteriori* dos parâmetros do VAR da volatilidade estocástica. A criação de funções separadas traz o ganho de se poder aplicar os modelos a grande quantidade de dados e salvar as saídas das simulações em objetos do tipo `ff`, que utilizam pouca memória RAM e, consequentemente, permitem armazenar grande quantidade de simulações.
  - Criar funções que amostram parâmetros de um modelo de regressão. A ideia é aproveitar a função que simula as cargas fatoriais, simplificando-a para o resultado dos modelos de regressão (similar ao caso s = 0).
  - Aplicar modelo de regressão fatorial e começar análise dos dados.
  - Transcrever funções do modelo dinâmico matriz-variado para o pacote.
  
#### 2014-06-21

Passos dados:

- Revistos os passos para a otimização da função de perda quadrática ex-post. Ainda não foi solucionado o problema. Pelos testes realizados, a otimização parece estar correta. Preciso ler o artigo com mais calma para tentar encontrar qual a nuância que deixei passar.

Próximos passos:

- Aplicar modelo de regressão fatorial para começar a analisar os dados;
- Rever códigos da otimização da função de perda;
- Juntar modelo de volatilidade estocástica aos modelos fatoriais na forma PLT.

#### 2014-06-19

Passos dados:

- Todas as funções para a aplicação do modelo fatorial dinâmico (`SampleDynFactors`, `SampleVarParms` e `SampleDfmLoads`) foram validadas, comparadas com a estimação via mínimos quadrados e avaliadas quanto à convergência da *posteriori* para o verdadeiro valor quando a amostra aumenta.
- A otimização da função de perda para solucionar o problema da rotação não foi solucionada. Avancei no estudo da otimização, mas os resultados ainda não estão corretos. Talvez haja algum erro na definição da função de perda ou alguma conta feita erroneamente.

Passos futuros:

- Revisar o artigo e encontrar a correta definição da função de perda quadrática que soluciona o problema de rotação.
- Aplicar modelo de regressão estática pelo Brent e por harmônicos onde a variância é definida por fatores latentes. Essa aplicação requer somente o uso de funções já definidas e a aplicação do WOP para solucionar o problema da identificabilidade da matriz de cargas.
- Aplicar modelo de regressão com volatilidade estocástica fatorial. Essa é uma extensão ao modelo anterior.

#### 2014-06-05

Passos dados:

- Funções `SampleDynFactors` e `SampleVarParms` criadas e validadas via simulação. Talvez ainda seja necessário rever o código da amostragem da distribuição *a posteriori* dos parâmetros do VAR. Parece que está um pouco frágil a estimação, mas pode ser uma dificuldade intrínseca do modelo.
- Abordagem WOP para matriz de cargas no DFM iniciada. Acabei não conseguindo terminar o código a tempo da reunião com a Alexandra ocorrida hoje. A questão é como otimizar a função de perda definida como soma de duas parcelas. A primeira parcela se refere ao WOP (bastante semelhante ao caso estático), enquanto a segunda se refere à perda quando aos parâmetros do VAR, que sofrem alteração devida à rotação dos fatores.

Passos futuros:

- Avançar no código para o WOP no caso dinâmico;
- Incluir SV aos modelos;
- Simulação das variâncias idiossincráticas seguindo SV sera feita utilizando a abordagem de Kastner e Frühwirth-Schnatter (2013), que é bastante eficiente para modelos univariados (caso em que cairemos);
- Construir funções que calculam a verossimilhança marginal para cada modelo. Ganha-se eficiência quanto mais integrais forem resolvidas analiticamente. Para isso, é interessante construir um filtro progressivo que possibilite integração dos parâmetros de estado e dos fatores estáticos.

#### 2014-05-20

Passos dados:

- Relatório de aplicação a dados simulados concluído. As funções estão funcionando corretamente e há convergência das cadeias, embora demore um pouco se o ponto de partida for muito distante do verdadeiro parâmetro. FFBS usando Choleski é mais rápido.

Passos futuros:

- Escrever relatório de aplicação a dados reais e observar comportamento das cadeias. Fazer comparação de modelos utilizando fator de Bayes;
- Incoporar SV aos modelos dinâmicos;
- Programar função para modelo fatorial dinâmico e DFSVM (modelo fatorial dinâmico de volatilidade estocástica).

#### 2014-05-19

Passos dados:

- Algoritmo de Gibbs implementado utilizando fatoração de Choleski para a geração dos vetores aleatórios normais no FFBS e também para a inversão de matrizes. Com essa abordagem houve ganho de 25% no tempo computacional para dados artificiais com $T = 500$ e $q = 9$.

Passos futuros:

- Concluir relatório de aplicação a dados artificiais;
- Incorporar SV aos modelos dinâmicos fatoriais;
- Programar modelo fatorial dinâmico.

#### 2014-05-18

Passos dados:

 - Cálculo da verosimilhança marginal, tendo integrado analiticamente tanto os fatores quanto os parâmetros de estado. Ponto a observar: a verossimilhança utilizando fatoração de Choleski antes da inversão parece ser mais eficiente do que a decomposição espectral (autovalores).
 
#### 2014-05-15

Passos dados:

  - tentativa de acelerar o MCMC via remocao do loop. Resultado: nao houve ganho computacional. Ao contrario, aumentou-se o tempo computacional.

Passos futuros:

  - implementar função de cálculo da verossimilhança marginal. Para tanto, é recomendável integrar tanto os fatores quanto os parâmetros de estado do MCMC;
  - aplicar o modelo a dados artificiais e gerar relatório dessa aplicação, ilustrando as diferenças entre os valores verdadeiros e os estimados pelo modelo. Atentar para a relação entre os parâmetros fixos e os dinâmicos;
  - incorporar SV aos modelos já implementados;
  - estudar modelos fatoriais dinâmicos (*dynamic factor models*), pois me parece uma abordagem interessante para comparar a análise fatorial com o modelo de regressão pelo Brent. Nessa abordagem, reduzir-se-ia a complexidade do modelo por não haver mais parâmetros dinâmicos na relação com o Brent e os níveis médio e sazonal seriam tratados constantes por simplificação. O ganho dessa abordagem estaria na interpretação que se poderia ter dos fatores globais, regionais e dos produtos.
