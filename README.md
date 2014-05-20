README
========================================================

***2014-05-19***

Passos dados:
- Algoritmo de Gibbs implementado utilizando fatoração de Choleski para a geração dos vetores aleatórios normais no FFBS e também para a inversão de matrizes. Com essa abordagem houve ganho de 25% no tempo computacional para dados artificiais com $T = 500$ e $q = 9$.

Passos futuros:
- Concluir relatório de aplicação a dados artificiais;
- Incorporar SV aos modelos dinâmicos fatoriais;
- Programar modelo fatorial dinâmico.

***2014-05-18***

Passos dados:
 - Cálculo da verosimilhança marginal, tendo integrado analiticamente tanto os fatores quanto os parâmetros de estado. Ponto a observar: a verossimilhança utilizando fatoração de Choleski antes da inversão parece ser mais eficiente do que a decomposição espectral (autovalores).
 
***2014-05-15***

Passos dados:
  - tentativa de acelerar o MCMC via remocao do loop. Resultado: nao houve ganho computacional. Ao contrario, aumentou-se o tempo computacional.

Passos futuros:
  - implementar função de cálculo da verossimilhança marginal. Para tanto, é recomendável integrar tanto os fatores quanto os parâmetros de estado do MCMC;
  - aplicar o modelo a dados artificiais e gerar relatório dessa aplicação, ilustrando as diferenças entre os valores verdadeiros e os estimados pelo modelo. Atentar para a relação entre os parâmetros fixos e os dinâmicos;
  - incorporar SV aos modelos já implementados;
  - estudar modelos fatoriais dinâmicos (*dynamic factor models*), pois me parece uma abordagem interessante para comparar a análise fatorial com o modelo de regressão pelo Brent. Nessa abordagem, reduzir-se-ia a complexidade do modelo por não haver mais parâmetros dinâmicos na relação com o Brent e os níveis médio e sazonal seriam tratados constantes por simplificação. O ganho dessa abordagem estaria na interpretação que se poderia ter dos fatores globais, regionais e dos produtos.
