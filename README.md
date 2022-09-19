# PF-ControladoresModuloLunar

Bem vindo ao repositório PF-ControladoresModuloLunar! Esse conjunto de códigos refere-se ao projeto final do curso de Engenharia de Controle e Automação, cursado pela aluna Laís Paulo Lima, que abordou o problema de aterrissagem lunar e cujo tema envolve o projeto de controladores para o sistema dinâmico de módulo lunar. 

Quatro scripts foram desenvolvidos para o controlador LQR: LM _PureDynamics.m, LM _DynamicsOptimized.m, LM _LQR1.m e LM _LQR2.m. O arquivo LM _PureDynamics.m deve ser rodado primeiro, pois nesse script é feita a representação do sistema dinâmico. Em segunda instância, deve-se rodar o script LM_DynamicsOptimized.m, onde constam a linearização do sistema e algumas análises importantes, como a observabilidade e a controlabilidade. Em seguida, os scripts seguintes podem ser rodados independentemente, uma vez que nos códigos LM _LQR1.m e LM _LQR2.m a construção do controlador LQR é feita em duas abordagens diferentes.

Três scripts foram desenvolvidos para o controlador LQG: LM _PureDynamicsLQG.m, LM _DynamicsOptimizedLQG.m e o LQG.m. A ordem de simulação deve ser exatamente esta. No primeiro código, a dinâmica é simulada; no segundo código, alguns cálculos importantes são feitos e, por fim, no script LQG.m o controlador LQG é projetado.
