#pragma once

#include<vector>
#include <string>
#include <iostream>
#include <random>

class GA
{
private:
	int max_genom_list;//個体数
	int item_num;//品物の個数
	int max_weight;//最大重量
	double crossoverRate = 0.8;//交叉率
	double individualMutationRate = 0.3;//個体突然変異率
	double genomMutaionRate = 0.3;//遺伝子突然変異率
public:
	std::vector<double> weight;//重さの配列
	std::vector<double> value;//価値の配列
	double resultSumValue;//評価関数の合計

	class Data//データ格納用クラス
	{
	private:
		double coefficient = 0.001;//評価関数用の定数
		int item_num;//品物の数
	public:
		std::vector<bool> isIncluded;//品物を入れるかどうかの配列
		double result;//評価関数の値
		double r_weight;//その個体の合計重さ
		double r_value;//その個体の合計価値

		Data(int _item_num)//コンストラクタ
		{
			item_num = _item_num;

			isIncluded.resize(item_num);//isIncludedの配列の長さの設定
		}

		double calcResult(double _maxWeight)//評価関数
		{
			if (r_weight > _maxWeight)//最大重量を超えている時
			{
				if (r_weight != 0)
					result = coefficient / ((r_weight-_maxWeight));//超えるほど評価関数の値が悪くなるように
				else
					result = 0;
				//result = r_value * coefficient;//これだと超えるほど評価関数の値が良くなる
			}
			else//最大重量を超えなければそのまま
			{
				result = r_value;
			}
			return result;
		}
	};

	std::vector<Data> data, prev_data;//操作前後で値を保持するために2個
	GA(int _max_genom_list, int _item_num, int _max_weight, std::vector<double> weight, std::vector<double> value);	//コンストラクタ
	bool init();//初期化
	bool selection();//選択
	bool crossover();//交叉
	bool mutation();//突然変異
	bool calc(bool enableDisplay);//評価関数の計算

	~GA();//デコンストラクタ
};


