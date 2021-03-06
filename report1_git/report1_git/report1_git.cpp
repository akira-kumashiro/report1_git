// report1.cpp : アプリケーションのエントリ ポイントを定義します。
//

#include "stdafx.h"
#include "GA.h"

#define MAX_GENERATION 30000
#define MAX_GENOM_LIST 50
#define MAX_WEIGHT 60
#define ITEM_NUM 50

int main()
{
	const double weight[ITEM_NUM] =  //品物の重さ
	{
		9, //1
		7, //2
		8, //3
		2, //4
		10, //5
		7, //6
		7, //7
		8, //8
		5, //9
		4, //10
		7, //11
		5, //12
		7, //13
		5, //14
		9, //15
		9, //16
		9, //17
		8, //18
		8, //19
		2, //20
		7, //21
		7, //22
		9, //23
		8, //24
		4, //25
		7, //26
		3, //27
		9, //28
		7,   //29
		7,   //30
		9, //31
		5, //32
		10, //33
		7, //34
		10, //35
		10, //36
		7, //37
		10, //38
		10, //39
		10, //40
		3, //41
		8, //42
		3, //43
		4, //44
		2, //45
		2, //46
		5, //47
		3, //48
		9,   //49
		2   //50
	};
	const double prise[ITEM_NUM] = //品物の値段
	{
		20, //1
		28, //2
		2, //3
		28, //4
		15, //5
		28, //6
		21, //7
		7, //8
		28, //9
		12, //10
		21, //11
		4, //12
		31, //13
		28, //14
		24, //15
		36, //16
		33, //17
		2, //18
		25, //19
		21, //20
		35, //21
		14, //22
		36, //23
		25, //24
		12, //25
		14, //26
		40, //27
		36, //28
		2,   //29
		28,   //30
		33, //31
		40, //32
		22, //33
		2, //34
		18, //35
		22, //36
		14, //37
		22, //38
		15, //39
		22, //40
		40, //41
		7, //42
		4, //43
		21, //44
		21, //45
		28, //46
		40, //47
		4, //48
		24,   //49
		21   //50
	};

	//配列をstd::vectorへ変換
	std::vector<double> w(weight, std::end(weight));
	std::vector<double> p(prise, std::end(prise));
	GA ga(MAX_GENOM_LIST, ITEM_NUM, MAX_WEIGHT, w, p);//遺伝的アルゴリズム諸関数をまとめたクラスの宣言

	ga.init();//変数の初期化

	for (int i = 0; i <= MAX_GENERATION; i++)//メインのループ
	{
		bool change = ga.selection();//選択
		ga.uniformityCrossover();//交叉

		ga.mutation();//突然変異
		
		if (i % (MAX_GENERATION / 10) == 0 || change)
		{
			std::cout << "i=" << std::to_string(i) << std::endl;
			ga.calc(true);//評価関数の計算
		}
	}
	std::sort(ga.data.begin(), ga.data.end(), [](const GA::Data& x, const GA::Data& y) { return x.result > y.result; });

	while (!_kbhit())
	{

	}
	return 0;
}

