#include "stdafx.h"
#include "GA.h"

GA::GA(int _max_genom_list, int _item_num, int _max_weight, std::vector<double> _weight, std::vector<double> _value) :
	data(std::vector<Data>(_max_genom_list, _item_num)),//data�̏�����
	eliteData(_item_num)
{
	//��������ϐ����N���X���ϐ��Ɋi�[
	max_genom_list = _max_genom_list;
	item_num = _item_num;
	max_weight = _max_weight;
	weight = _weight;
	value = _value;

	//�z��̒����̐ݒ�
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
			data[i].isIncluded[j] = (random(0,1) == 1 ? true : false);//��`�q�̏����ݒ�
			//data[i].isIncluded[j] = false;//�����l�S��0�̂Ƃ�
			printf_s("%d", data[i].isIncluded[j] ? 1 : 0);

			data[i].r_value += value[j] * data[i].isIncluded[j];
			data[i].r_weight += weight[j] * data[i].isIncluded[j];
			//�d���ƒl�i�̌v�Z
		}
		data[i].calcResult(max_weight);//�]���֐�
		resultSumValue += data[i].result;//�]���֐��̍��v���v�Z
		printf_s(" \t sumValue=%4.0lf\t sumWeight=%4.0lf\t Result=%7.4lf\n", data[i].r_value, data[i].r_weight, data[i].result);
	}
	return true;
}

bool GA::selection()
{
	int max_num = 0;//�ł��]���̗ǂ��̂̔ԍ�
	bool ret = false;

	calc(false);

	prev_data = data;

	resultSumValue = 0;
	for (int i = 0; i < max_genom_list; i++)
		//���[���b�g�I��p�ɕ]���֐��̍��v�ƈ�ԕ]���̗ǂ��ԍ����擾
	{
		resultSumValue += prev_data[i].result;
		if (prev_data[i].result > prev_data[max_num].result)
			max_num = i;
	}

	eliteData = prev_data[max_num];//�f�[�^�̐擪�͍ł��]���̗ǂ���
//	eliteData = prev_data[minNum];
	if (eliteData.result - prev_data[minNum].result >= 1)//�ł��]���̗ǂ��̂̕ω��̊Ď�(�f�o�b�O�p)
		ret = true;

	for (int i = 0; i < max_genom_list; i++)
	{
		double selector = random(0.0,1.0);//�����𐶐�
		double needle = 0;//���[���b�g�̐j�𐶐�
		int j = 0;
		for (;; j++)
		{
			needle += (prev_data[j].result / resultSumValue);//���[���b�g�̐j�𗐐��̒l�܂Ői�߂�
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

	for (int i = 0; i < max_genom_list; i += 2)//2������
	{
		for (int j = 0; j < item_num; j++)
		{
			bool isCrossover = (random(0.0,1.0) >= crossoverRate ? true : false);//true�Ō����Ȃ�
			data[i + 1].isIncluded[j] = isCrossover ? prev_data[i + 1].isIncluded[j] : prev_data[i].isIncluded[j];
			//if (i != 0)//�擪�̃f�[�^�͕ی�
			data[i].isIncluded[j] = isCrossover ? prev_data[i].isIncluded[j] : prev_data[i + 1].isIncluded[j];
		}
	}
	return true;
}

bool GA::onePointCrossover()
{
	prev_data = data;

	for (int i = 0; i < max_genom_list; i += 2)//2������
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

	for (int i = 0; i < max_genom_list; i += 2)//2������
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

	for (int i = 0; i < max_genom_list; i += 2)//2������
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
		if (random(0.0, 1.0) <= individualMutationRate)//�̓ˑR�ψٗ��̌v�Z
		{
			int point = random(0, item_num - 1);
			data[i].isIncluded[point] = !data[i].isIncluded[point];
			//			for (int j = 0; j < item_num; j++)
			//			{
			//				data[i].isIncluded[j] = distribution(engine) <= genomMutaionRate ? !data[i].isIncluded[j] : data[i].isIncluded[j];//��`�q�ˑR�ψٗ��̌v�Z�@�ω�����ꍇ�̓r�b�g���]
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
			data[i].r_value += value[j] * data[i].isIncluded[j];//�̂̍��v���l���v�Z
			data[i].r_weight += weight[j] * data[i].isIncluded[j];//�̂̍��v�d�����v�Z
		}
		data[i].calcResult(max_weight);//�]���֐��̌v�Z
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
				printf_s("%d", data[i].isIncluded[j] ? 1 : 0);//�f�o�b�O�p
			}
			printf_s(" \t sumValue=%.0lf\t sumWeight=%.0lf\t Result=%.4lf\n", data[i].r_value, data[i].r_weight, data[i].result);
		}
	}

	if (enableDisplay)//�f�o�b�O�p
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
	//�����̐ݒ�
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
