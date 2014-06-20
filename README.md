README
========================================================

#### 2014-06-19

Passos dados:

- Todas as funções para a aplicação do modelo fatorial dinâmico (`SampleDynFactors`, `SampleVarParms` e `SampleDfmLoads`) foram validadas, comparadas com a estimação via mínimos quadrados e avaliadas quanto a convergência da posteriori para o verdadeiro valor quando a amostra aumenta.
- A otimização da função de perda para solucionar o problema da rotação não foi solucionada. Avancei no estudo da otimização, mas os resultados ainda não estão corretos. Talvez haja algum erro de definição da função de perda ou alguma conta feita erroneamente.

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
