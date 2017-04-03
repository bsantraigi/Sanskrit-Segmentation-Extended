import pickle
import json
class MatDB:
	def __init__(self):
		self.mat_lemCount_1D = pickle.load(open('../NewData/gauravs/Temporary_1D/mat_lemCount_1D.p', 'rb'), encoding = u'utf-8')
		self.mat_cngCount_1D = pickle.load(open('../NewData/gauravs/Temporary_1D/mat_cngCount_1D.p', 'rb'), encoding = u'utf-8')
		self.mat_tupCount_1D = pickle.load(open('../NewData/gauravs/Temporary_1D/mat_tupCount_1D.p', 'rb'), encoding = u'utf-8')

		self.mat_lemCount_1D = pickle.load(open('../NewData/gauravs/1D/mat_lemCount_1D.p', 'rb'), encoding = u'utf-8')
		self.mat_cngCount_1D = pickle.load(open('../NewData/gauravs/1D/mat_cngCount_1D.p', 'rb'), encoding = u'utf-8')
		self.mat_tupCount_1D = json.load(open('../NewData/gauravs/1D/mat_tupCount_1D.json'))

		self.mat_lem2lem_countonly = pickle.load(open('../NewData/gauravs/mat_lem2lem_old_countonly.p', 'rb'), encoding = u'utf-8')
		self.mat_lem2cng_countonly = pickle.load(open('../NewData/gauravs/mat_lemma2cng_new_countonly.p', 'rb'), encoding = u'utf-8')
		self.mat_lem2tup_countonly = pickle.load(open('../NewData/gauravs/mat_lem2tup_old_countonly.p', 'rb'), encoding = u'utf-8')

		self.mat_cng2lem_countonly = pickle.load(open('../NewData/gauravs/mat_cng2lemma_new_countonly.p', 'rb'), encoding = u'utf-8')
		self.mat_cng2tup_countonly = pickle.load(open('../NewData/gauravs/mat_cng2tup_new_countonly.p', 'rb'), encoding = u'utf-8')
		self.mat_cng2cng_countonly = pickle.load(open('../NewData/gauravs/mat_cng2cng_new_countonly.p', 'rb'), encoding = u'utf-8')
		
		self.mat_tup2cng_countonly = pickle.load(open('../NewData/gauravs/mat_tup2cng_new_countonly.p', 'rb'), encoding = u'utf-8')
		self.mat_tup2lem_countonly = pickle.load(open('../NewData/gauravs/mat_tup2lem_old_countonly.p', 'rb'), encoding = u'utf-8')
		self.mat_tup2tup_countonly = pickle.load(open('../NewData/gauravs/mat_tup2tup_new_countonly.p', 'rb'), encoding = u'utf-8')

		# self.prob_lem2lem_linear = pickle.load(open('../NewData/gauravs/probmat_linear/prob_lem2lem_linear.p', 'rb'), encoding = u'utf-8')
		# self.prob_lem2cng_linear = pickle.load(open('../NewData/gauravs/probmat_linear/prob_lem2cng_linear.p', 'rb'), encoding = u'utf-8')
		# self.prob_lem2tup_linear = pickle.load(open('../NewData/gauravs/probmat_linear/prob_lem2tup_linear.p', 'rb'), encoding = u'utf-8')

		# self.prob_cng2lem_linear = pickle.load(open('../NewData/gauravs/probmat_linear/prob_cng2lem_linear.p', 'rb'), encoding = u'utf-8')
		# self.prob_cng2tup_linear = pickle.load(open('../NewData/gauravs/probmat_linear/prob_cng2tup_linear.p', 'rb'), encoding = u'utf-8')
		# self.prob_cng2cng_linear = pickle.load(open('../NewData/gauravs/probmat_linear/prob_cng2cng_linear.p', 'rb'), encoding = u'utf-8')

		# self.prob_tup2cng_linear = pickle.load(open('../NewData/gauravs/probmat_linear/prob_tup2cng_linear.p', 'rb'), encoding = u'utf-8')
		# self.prob_tup2lem_linear = pickle.load(open('../NewData/gauravs/probmat_linear/prob_tup2lem_linear.p', 'rb'), encoding = u'utf-8')
		# self.prob_tup2tup_linear = pickle.load(open('../NewData/gauravs/probmat_linear/prob_tup2tup_linear.p', 'rb'), encoding = u'utf-8')
