#pragma once

#include<vector>
#include <string>
#include <iostream>
#include <random>

class GA
{
private:
	int max_genom_list;//�̐�
	int item_num;//�i���̌�
	int max_weight;//�ő�d��
	double crossoverRate = 0.8;//������
	double individualMutationRate = 0.3;//�̓ˑR�ψٗ�
	double genomMutaionRate = 0.3;//��`�q�ˑR�ψٗ�
public:
	std::vector<double> weight;//�d���̔z��
	std::vector<double> value;//���l�̔z��
	double resultSumValue;//�]���֐��̍��v

	class Data//�f�[�^�i�[�p�N���X
	{
	private:
		double coefficient = 0.001;//�]���֐��p�̒萔
		int item_num;//�i���̐�
	public:
		std::vector<bool> isIncluded;//�i�������邩�ǂ����̔z��
		double result;//�]���֐��̒l
		double r_weight;//���̌̂̍��v�d��
		double r_value;//���̌̂̍��v���l

		Data(int _item_num)//�R���X�g���N�^
		{
			item_num = _item_num;

			isIncluded.resize(item_num);//isIncluded�̔z��̒����̐ݒ�
		}

		double calcResult(double _maxWeight)//�]���֐�
		{
			if (r_weight > _maxWeight)//�ő�d�ʂ𒴂��Ă��鎞
			{
				if (r_weight != 0)
					result = coefficient / ((r_weight-_maxWeight));//������قǕ]���֐��̒l�������Ȃ�悤��
				else
					result = 0;
				//result = r_value * coefficient;//���ꂾ�ƒ�����قǕ]���֐��̒l���ǂ��Ȃ�
			}
			else//�ő�d�ʂ𒴂��Ȃ���΂��̂܂�
			{
				result = r_value;
			}
			return result;
		}
	};

	std::vector<Data> data, prev_data;//����O��Œl��ێ����邽�߂�2��
	GA(int _max_genom_list, int _item_num, int _max_weight, std::vector<double> weight, std::vector<double> value);	//�R���X�g���N�^
	bool init();//������
	bool selection();//�I��
	bool crossover();//����
	bool mutation();//�ˑR�ψ�
	bool calc(bool enableDisplay);//�]���֐��̌v�Z

	~GA();//�f�R���X�g���N�^
};


