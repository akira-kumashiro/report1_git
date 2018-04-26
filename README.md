# report1_git
アルゴリズム特論課題1(C++)  
C++で書きました
## 概略
`GA`クラスに遺伝的アルゴリズム関係の関数をまとめました  
`GA::Data`に遺伝子と評価値を格納してます

## 交叉の種類
### 1. [`uniformityCrossover()`](https://github.com/akira-kumashiro/report1_git/blob/master/report1_git/report1_git/GA.cpp#L84 "uniformityCrossover()")
一様交叉  
乱数生成が多くて遅い
### 2. [`onePointCrossover()`](https://github.com/akira-kumashiro/report1_git/blob/master/report1_git/report1_git/GA.cpp#L101 "onePointCrossover()")
0~49の中でどれか1点を決めてそこを境にブロック分け  
ブロックごとに交叉率分交叉する
### 3. [`twoPointCrossover()`](https://github.com/akira-kumashiro/report1_git/blob/master/report1_git/report1_git/GA.cpp#L128 "twoPointCrossover()")
`onePointCrossover()`のブロックが3つになった版
### 4. [`tsunoPointCrossover()`](https://github.com/akira-kumashiro/report1_git/blob/master/report1_git/report1_git/GA.cpp#L165 "tsunoPointCrossover()")
2点交叉だけど点を打つのは0~49の2点  
2点の間で交叉を行う  

## 突然変異について
突然変異を起こしたとき書き換えられるのは1点だけ
