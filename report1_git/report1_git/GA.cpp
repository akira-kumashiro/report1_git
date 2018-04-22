#include "stdafx.h"
#include "GA.h"

GA::GA(int _max_genom_list, int _item_num, int _max_weight, std::vector<double> _weight, std::vector<double> _value) :
	data(std::vector<Data>(_max_genom_list, _item_num))//data‚Ì‰Šú‰»
{
	//‚à‚ç‚Á‚½•Ï”‚ğƒNƒ‰ƒX“à•Ï”‚ÉŠi”[
	max_genom_list = _max_genom_list;
	item_num = _item_num;
	max_weight = _max_weight;
	weight = _weight;
	value = _value;

	//”z—ñ‚Ì’·‚³‚Ìİ’è
	weight.resize(item_num);
	value.resize(item_num);
}

bool GA::init()
{
	//—”‚Ìİ’è
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_int_distribution<int> distribution(0, 1);

	for (int i = 0; i < max_genom_list; i++)
	{
		data[i].r_value = 0;
		data[i].r_weight = 0;
		for (int j = 0; j < item_num; j++)
		{
			data[i].isIncluded[j] = (distribution(engine) == 1 ? true : false);//ˆâ“`q‚Ì‰Šúİ’è
																			   //data[i].isIncluded[j] = false;
			printf_s("%d", data[i].isIncluded[j] ? 1 : 0);

			data[i].r_value += value[j] * data[i].isIncluded[j];
			data[i].r_weight += weight[j] * data[i].isIncluded[j];
			//d‚³‚Æ’l’i‚ÌŒvZ
		}
		data[i].calcResult(max_weight);//•]‰¿ŠÖ”
		resultSumValue += data[i].result;//•]‰¿ŠÖ”‚Ì‡Œv‚ğŒvZ
		printf_s(" \t sumValue=%4.0lf\t sumWeight=%4.0lf\t Result=%7.4lf\n", data[i].r_value, data[i].r_weight, data[i].result);
	}
	return true;
}

bool GA::selection()
{
	int max_num = 0;//Å‚à•]‰¿‚Ì—Ç‚¢ŒÂ‘Ì‚Ì”Ô†
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_real_distribution<double> distribution(0, 1.0);
	bool ret = false;
	prev_data = data;

	resultSumValue = 0;
	for (int i = 0; i < max_genom_list; i++)
	{
		resultSumValue += prev_data[i].result;
		if (prev_data[i].result > prev_data[max_num].result)
		{
			max_num = i;
		}
	}

	data[0] = prev_data[max_num];
	if (data[0].result - prev_data[0].result >= 1)
	{
		ret = true;
	}

	for (int i = 1; i < max_genom_list; i++)
	{
		double selector = distribution(engine);
		double needle = 0;
		int j = 0;
		for (;; j++)
		{
			needle += prev_data[j].result / resultSumValue;
			if (needle > selector)
				break;
			if (j == max_genom_list - 1)
				break;
		}
		data[i] = prev_data[j];
	}

	return ret;
}

bool GA::crossover()
{
	prev_data = data;

	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_real_distribution<double> distribution(0, 1.0);
	std::vector<bool> uniform_temp(item_num);

	for (int i = 0; i < max_genom_list; i += 2)
	{
		for (int j = 0; j < item_num; j++)
		{
			uniform_temp[j] = (distribution(engine) >= crossoverRate ? true : false);//true‚ÅŒğ³‚È‚µ
		}
		for (int j = 0; j < item_num; j++)
		{
			data[i + 1].isIncluded[j] = uniform_temp[j] ? prev_data[i + 1].isIncluded[j] : prev_data[i].isIncluded[j];
			if (j == 0)
				break;
			data[i].isIncluded[j] = uniform_temp[j] ? prev_data[i].isIncluded[j] : prev_data[i + 1].isIncluded[j];
		}
	}

	return true;
}

bool GA::mutation()
{
	std::mt19937 engine;
	std::uniform_real_distribution<double> distribution(0, 1.0);
	std::vector<bool> uniform_temp(item_num);

	for (int i = 1; i < max_genom_list; i++)
	{
		if (distribution(engine) <= individualMutationRate)
		{
			for (int j = 0; j < item_num; j++)
			{
				data[i].isIncluded[j] = distribution(engine) <= genomMutaionRate ? !data[i].isIncluded[j] : data[i].isIncluded[j];
			}
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
			data[i].r_value += value[j] * data[i].isIncluded[j];
			data[i].r_weight += weight[j] * data[i].isIncluded[j];

			if (enableDisplay)
			{
				printf_s("%d", data[i].isIncluded[j] ? 1 : 0);
			}
		}
		data[i].calcResult(max_weight);
		if (enableDisplay)
			printf_s(" \t sumValue=%.0lf\t sumWeight=%.0lf\t Result=%.4lf\n", data[i].r_value, data[i].r_weight, data[i].result);
	}

	if (enableDisplay)
	{
		for (int i = 0; i < item_num; i++)
		{
			if (data[0].isIncluded[i])
				printf_s("+%d", (int)value[i]);
		}
		printf_s("=%lf\n", data[0].result);
	}

	return true;
}

GA::~GA()
{

}
