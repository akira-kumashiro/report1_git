#include "stdafx.h"
#include "GA.h"

GA::GA(int _max_genom_list, int _item_num, int _max_weight, std::vector<double> _weight, std::vector<double> _value) :
	data(std::vector<Data>(_max_genom_list, _item_num))//data�̏�����
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
	//�����̐ݒ�
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_int_distribution<int> distribution(0, 1);

	for (int i = 0; i < max_genom_list; i++)
	{
		data[i].r_value = 0;
		data[i].r_weight = 0;
		for (int j = 0; j < item_num; j++)
		{
			data[i].isIncluded[j] = (distribution(engine) == 1 ? true : false);//��`�q�̏����ݒ�
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
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_real_distribution<double> distribution(0, 1.0);
	bool ret = false;
	prev_data = data;

	resultSumValue = 0;
	for (int i = 0; i < max_genom_list; i++)
		//���[���b�g�I��p�ɕ]���֐��̍��v�ƈ�ԕ]���̗ǂ��ԍ����擾
	{
		resultSumValue += prev_data[i].result;
		if (prev_data[i].result > prev_data[max_num].result)
			max_num = i;
	}

	data[0] = prev_data[max_num];//�f�[�^�̐擪�͍ł��]���̗ǂ���
	if (data[0].result - prev_data[0].result >= 1)//�ł��]���̗ǂ��̂̕ω��̊Ď�(�f�o�b�O�p)
		ret = true;

	for (int i = 1; i < max_genom_list; i++)
	{
		double selector = distribution(engine);//�����𐶐�
		double needle = 0;//���[���b�g�̐j�𐶐�
		int j = 0;
		for (;; j++)
		{
			needle += prev_data[j].result / resultSumValue;//���[���b�g�̐j�𗐐��̒l�܂Ői�߂�
			if (needle > selector)
				break;
			if (j == max_genom_list - 1)
				break;
		}
		data[i] = prev_data[j];
	}
	return ret;
}

bool GA::uniformityCrossover()
{
	prev_data = data;

	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_real_distribution<double> distribution(0, 1.0);

	for (int i = 0; i < max_genom_list; i += 2)//2������
	{
		for (int j = 0; j < item_num; j++)
		{
			bool isCrossover = (distribution(engine) >= crossoverRate ? true : false);//true�Ō����Ȃ�
			data[i + 1].isIncluded[j] = isCrossover ? prev_data[i + 1].isIncluded[j] : prev_data[i].isIncluded[j];
			if (i == 0)//�擪�̃f�[�^�͕ی�
				break;
			data[i].isIncluded[j] = isCrossover ? prev_data[i].isIncluded[j] : prev_data[i + 1].isIncluded[j];
		}
	}
	return true;
}

bool GA::onePointCrossover()
{
	prev_data = data;

	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_real_distribution<double> distribution(0, 1.0);

	for (int i = 0; i < max_genom_list; i += 2)//2������
	{
		if (distribution(engine) <= crossoverRate)
		{
			for (int j = item_num / 2 - 1; j < item_num; j++)
			{
				//bool isCrossover = (distribution(engine) >= crossoverRate ? true : false);//true�Ō����Ȃ�
				data[i + 1].isIncluded[j] = prev_data[i].isIncluded[j];
				if (i == 0)//�擪�̃f�[�^�͕ی�
					break;
				data[i].isIncluded[j] = prev_data[i + 1].isIncluded[j];
			}
		}
	}
	return true;
}

bool GA::twoPointCrossover()
{
	prev_data = data;

	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_real_distribution<double> distribution(0, 1.0);
	//std::uniform_int_distribution<int> disInt(0, 25);


	for (int i = 0; i < max_genom_list; i += 2)//2������
	{
		/*double d1 = distribution(engine) / 2;
		int del1 = d1*item_num;
		int del2 = d1 != 0 ? (distribution(engine) / (1 - d1) + d1)*item_num:distribution(engine)*item_num;*/
		std::uniform_int_distribution<int> disInt1(0, item_num / 2);
		int del1 = disInt1(engine);
		std::uniform_int_distribution<int> disInt2(del1, item_num - 1);
		int del2 = disInt2(engine);
		if (distribution(engine) <= crossoverRate)
		{
			for (int j = 0; j < del1; j++)
			{
				data[i + 1].isIncluded[j] = prev_data[i].isIncluded[j];
				if (i == 0)//�擪�̃f�[�^�͕ی�
					break;
				data[i].isIncluded[j] = prev_data[i + 1].isIncluded[j];
			}
		}
		if (distribution(engine) <= crossoverRate)
		{
			for (int j = del1; j < del2; j++)
			{
				data[i + 1].isIncluded[j] = prev_data[i].isIncluded[j];
				if (i == 0)//�擪�̃f�[�^�͕ی�
					break;
				data[i].isIncluded[j] = prev_data[i + 1].isIncluded[j];
			}
		}

		if (distribution(engine) <= crossoverRate)
		{
			for (int j = del2; j < item_num; j++)
			{
				data[i + 1].isIncluded[j] = prev_data[i].isIncluded[j];
				if (i == 0)//�擪�̃f�[�^�͕ی�
					break;
				data[i].isIncluded[j] = prev_data[i + 1].isIncluded[j];
			}
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
		if (distribution(engine) <= individualMutationRate)//�̓ˑR�ψٗ��̌v�Z
		{
			for (int j = 0; j < item_num; j++)
			{
				data[i].isIncluded[j] = distribution(engine) <= genomMutaionRate ? !data[i].isIncluded[j] : data[i].isIncluded[j];//��`�q�ˑR�ψٗ��̌v�Z�@�ω�����ꍇ�̓r�b�g���]
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
			data[i].r_value += value[j] * data[i].isIncluded[j];//�̂̍��v���l���v�Z
			data[i].r_weight += weight[j] * data[i].isIncluded[j];//�̂̍��v�d�����v�Z

			if (enableDisplay)
				printf_s("%d", data[i].isIncluded[j] ? 1 : 0);//�f�o�b�O�p
		}
		data[i].calcResult(max_weight);//�]���֐��̌v�Z
		if (enableDisplay)//�f�o�b�O�p
			printf_s(" \t sumValue=%.0lf\t sumWeight=%.0lf\t Result=%.4lf\n", data[i].r_value, data[i].r_weight, data[i].result);
	}

	if (enableDisplay)//�f�o�b�O�p
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
