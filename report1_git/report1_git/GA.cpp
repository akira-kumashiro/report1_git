#include "stdafx.h"
#include "GA.h"

GA::GA(int _max_genom_list, int _item_num, int _max_weight, std::vector<double> _weight, std::vector<double> _value) :
	data(std::vector<Data>(_max_genom_list, _item_num)),//dataの初期化
	eliteData(_item_num)
{
	//もらった変数をクラス内変数に格納
	max_genom_list = _max_genom_list;
	item_num = _item_num;
	max_weight = _max_weight;
	weight = _weight;
	value = _value;

	//配列の長さの設定
	weight.resize(item_num);
	value.resize(item_num);
}

bool GA::init()
{
	for (int i = 0; i < max_genom_list; i++)
	{
		data[i].r_value = 0;
		data[i].r_weight = 0;
		for (int j = 0; j < item_num; j++)
		{
			data[i].isIncluded[j] = (random(0,1) == 1 ? true : false);//遺伝子の初期設定
			//data[i].isIncluded[j] = false;//初期値全部0のとき
			printf_s("%d", data[i].isIncluded[j] ? 1 : 0);

			data[i].r_value += value[j] * data[i].isIncluded[j];
			data[i].r_weight += weight[j] * data[i].isIncluded[j];
			//重さと値段の計算
		}
		data[i].calcResult(max_weight);//評価関数
		resultSumValue += data[i].result;//評価関数の合計を計算
		printf_s(" \t sumValue=%4.0lf\t sumWeight=%4.0lf\t Result=%7.4lf\n", data[i].r_value, data[i].r_weight, data[i].result);
	}
	return true;
}

bool GA::selection()
{
	int max_num = 0;//最も評価の良い個体の番号
	bool ret = false;

	calc(false);

	prev_data = data;

	resultSumValue = 0;
	for (int i = 0; i < max_genom_list; i++)
		//ルーレット選択用に評価関数の合計と一番評価の良い番号を取得
	{
		resultSumValue += prev_data[i].result;
		if (prev_data[i].result > prev_data[max_num].result)
			max_num = i;
	}

	eliteData = prev_data[max_num];//データの先頭は最も評価の良い個体
//	eliteData = prev_data[minNum];
	if (eliteData.result - prev_data[minNum].result >= 1)//最も評価の良い個体の変化の監視(デバッグ用)
		ret = true;

	for (int i = 0; i < max_genom_list; i++)
	{
		double selector = random(0.0,1.0);//乱数を生成
		double needle = 0;//ルーレットの針を生成
		int j = 0;
		for (;; j++)
		{
			needle += (prev_data[j].result / resultSumValue);//ルーレットの針を乱数の値まで進める
			if (needle > selector)
				break;
			if (j == (max_genom_list - 1))
				break;
		}
		data[i] = prev_data[j];
	}
	return ret;
}

bool GA::uniformityCrossover()
{
	prev_data = data;

	for (int i = 0; i < max_genom_list; i += 2)//2個ずつ交叉
	{
		for (int j = 0; j < item_num; j++)
		{
			bool isCrossover = (random(0.0,1.0) >= crossoverRate ? true : false);//trueで交叉なし
			data[i + 1].isIncluded[j] = isCrossover ? prev_data[i + 1].isIncluded[j] : prev_data[i].isIncluded[j];
			//if (i != 0)//先頭のデータは保護
			data[i].isIncluded[j] = isCrossover ? prev_data[i].isIncluded[j] : prev_data[i + 1].isIncluded[j];
		}
	}
	return true;
}

bool GA::onePointCrossover()
{
	prev_data = data;

	for (int i = 0; i < max_genom_list; i += 2)//2個ずつ交叉
	{
		int del1 = random(0, item_num - 1);
		if (random(0.0,1.0) <= crossoverRate)
		{
			for (int j = 0; j < del1; j++)
			{
				data[i + 1].isIncluded[j] = prev_data[i].isIncluded[j];
				data[i].isIncluded[j] = prev_data[i + 1].isIncluded[j];
			}
		}
		if (random(0.0,1.0) <= crossoverRate)
		{
			for (int j = del1; j < item_num; j++)
			{
				data[i + 1].isIncluded[j] = prev_data[i].isIncluded[j];
				data[i].isIncluded[j] = prev_data[i + 1].isIncluded[j];
			}
		}
	}
	return true;
}

bool GA::twoPointCrossover()
{
	prev_data = data;

	for (int i = 0; i < max_genom_list; i += 2)//2個ずつ交叉
	{
		int del1 = random(0, item_num / 2);
		int del2 = random(del1, item_num - 1);
		if (random(0.0,1.0) <= crossoverRate)
		{
			for (int j = 0; j < del1; j++)
			{
				data[i + 1].isIncluded[j] = prev_data[i].isIncluded[j];
				data[i].isIncluded[j] = prev_data[i + 1].isIncluded[j];
			}
		}
		if (random(0.0, 1.0) <= crossoverRate)
		{
			for (int j = del1; j < del2; j++)
			{
				data[i + 1].isIncluded[j] = prev_data[i].isIncluded[j];
				data[i].isIncluded[j] = prev_data[i + 1].isIncluded[j];
			}
		}

		if (random(0.0, 1.0) <= crossoverRate)
		{
			for (int j = del2; j < item_num; j++)
			{
				data[i + 1].isIncluded[j] = prev_data[i].isIncluded[j];
				data[i].isIncluded[j] = prev_data[i + 1].isIncluded[j];
			}
		}
	}
	return true;
}

bool GA::tsunoPointCrossover()
{
	prev_data = data;

	for (int i = 0; i < max_genom_list; i += 2)//2個ずつ交叉
	{
		int del1 = random(0, item_num - 1);
		int del2 = random(del1, item_num);
		if (random(0.0,1.0) <= crossoverRate)
		{
			for (int j = del1; j < del2; j++)
			{
				data[i + 1].isIncluded[j] = prev_data[i].isIncluded[j];
				data[i].isIncluded[j] = prev_data[i + 1].isIncluded[j];
			}
		}
	}

	return true;
}

bool GA::mutation()
{
	for (int i = 0; i < max_genom_list; i++)
	{
		if (random(0.0, 1.0) <= individualMutationRate)//個体突然変異率の計算
		{
			int point = random(0, item_num - 1);
			data[i].isIncluded[point] = !data[i].isIncluded[point];
			//			for (int j = 0; j < item_num; j++)
			//			{
			//				data[i].isIncluded[j] = distribution(engine) <= genomMutaionRate ? !data[i].isIncluded[j] : data[i].isIncluded[j];//遺伝子突然変異率の計算　変化する場合はビット反転
			//			}
		}
	}
	return true;
}

bool GA::calc(bool enableDisplay)
{
	for (int i = 0; i < max_genom_list; i++)
	{
		data[i].r_value = 0;
		data[i].r_weight = 0;

		for (int j = 0; j < item_num; j++)
		{
			data[i].r_value += value[j] * data[i].isIncluded[j];//個体の合計価値を計算
			data[i].r_weight += weight[j] * data[i].isIncluded[j];//個体の合計重さを計算
		}
		data[i].calcResult(max_weight);//評価関数の計算
		if (data[i].result < data[minNum].result)
		{
			minNum = i;
		}
		if (data[i].result > data[maxNum].result)
		{
			maxNum = i;
		}
	}
	data[minNum] = eliteData;

	if (enableDisplay)
	{
		for (int i = 0; i < max_genom_list; i++)
		{
			for (int j = 0; j < item_num; j++)
			{
				printf_s("%d", data[i].isIncluded[j] ? 1 : 0);//デバッグ用
			}
			printf_s(" \t sumValue=%.0lf\t sumWeight=%.0lf\t Result=%.4lf\n", data[i].r_value, data[i].r_weight, data[i].result);
		}
	}

	if (enableDisplay)//デバッグ用
	{
		for (int i = 0; i < item_num; i++)
		{
			if (data[minNum].isIncluded[i])
				printf_s("+%d", (int)value[i]);
		}
		printf_s("=%lf\n", data[minNum].result);
	}
	return true;
}

int GA::random(int min, int max)
{
	//乱数の設定
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_int_distribution<int> distribution(min, max);
	return distribution(engine);
}

double GA::random(int min, double max)
{
	return random((double)min,max);
}

double GA::random(double min, int max)
{
	return random(min,(double)max);
}

double GA::random(double min, double max)
{
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(engine);
}



GA::~GA()
{

}
