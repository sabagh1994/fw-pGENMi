import os
import torch
import random
import numpy as np
import torch.nn as nn

def get_file_dims(input_file, header=True, **kwargs):
    kwargs.setdefault('output_file', None)
    num_columns = 0
    num_genes = 0
    with open(input_file) as fp:
        line_no = 0
        for line in fp:
            line_no += 1
            if (line_no == 1 and header == True):
                next
            else:
                new_num_columns = len(line.split("\t"))
                if (num_columns != 0 and num_columns != new_num_columns):
                    print("Error the number of columns isn't consistent", file=kwargs['output_file'])
                    print(line, file=kwargs['output_file']);
                    print("Num genes " + str(num_genes), file=kwargs['output_file'])
                    print("New columns " + str(new_num_columns), file=kwargs['output_file'])
                    print("Old columns " + str(num_columns), file=kwargs['output_file'])
                    exit(1)
                num_columns = new_num_columns
                num_genes += 1 
    return num_genes, num_columns
  
def load_data(input_file, num_features):
    evidence = np.loadtxt(input_file, dtype="float64", skiprows = 1, delimiter="\t", usecols=tuple([i for i in range(2, 2 + num_features)]))
    data = np.loadtxt(input_file, dtype="float64", skiprows = 1, delimiter="\t", usecols=1)
    genes = np.loadtxt(input_file, dtype="str", skiprows = 1, delimiter="\t", usecols=0)
    return genes, data, evidence
    
class PGM_LL_FW(nn.Module):
    '''
    Model class for factorized-weighted pgenmi (fwpgenmi)
    There are separate weights for TFs and dynamic marks
    '''
    def __init__(self, tf_dim, mark_dim):
        super(PGM_LL_modif, self).__init__()
        self.pgm_w_tf = torch.nn.Parameter(data=torch.empty(tf_dim, 1, dtype=torch.double))
        self.pgm_w_m = torch.nn.Parameter(data=torch.empty(1, mark_dim, dtype=torch.double))
        self.pgm_bias = torch.nn.Parameter(data=torch.empty(1, 1, dtype=torch.double))
        self.alpha= torch.nn.Parameter(data=torch.empty(1, dtype=torch.double), requires_grad= True)
      
    def forward(self, x):
        y = x@self.pgm_w
        return y
    
    @property
    def pgm_w(self):
        pgm_w_ = (self.pgm_w_tf@self.pgm_w_m).reshape(-1,1)
        pgm_w_ = torch.cat((self.pgm_bias, pgm_w_), dim= 0)
        return pgm_w_
    
class PGM_LL(nn.Module):
    '''
    Model class for pgenmi
    There are (TF,dynamic mark)-specific weights
    '''
    def __init__(self, in_dim, out_dim):
        super(PGM_LL, self).__init__()
        self.pgm_w = torch.nn.Parameter(data=torch.empty(in_dim, out_dim, dtype=torch.double))
        self.alpha = torch.nn.Parameter(data=torch.empty(1, dtype= torch.double), requires_grad=True)
    def forward(self, x):
        y = x@self.pgm_w
        return y
    
def estep(model, x, pg, regul=0, is_H0=False):
    prob_zg1 = torch.sigmoid(torch.squeeze(model(x))) #17511 is the shape
    a= model.alpha
    # Compute p(pg, zg=1 | xg1..xgm, w1...wm)
    prob_pg1 = a * torch.pow(pg, a - 1) * prob_zg1
    # Compute p(pg, zg=0 | xg1..xgm, w1...wm)
    prob_pg0 = 1 - prob_zg1
    # Compute p(pg| xg1..xgm, w1...wm)
    prob_pg = prob_pg1 + prob_pg0
    # Compute p(zg=1|pg, xg1..xgm, w1...wm)
    post_zg1 = torch.div(prob_pg1, prob_pg)
    # Compute total loglik
    if not is_H0:
        L2_loss= torch.sum(torch.dot(model.pgm_w[1:,:].reshape(-1), model.pgm_w[1:,:].reshape(-1)))
    else:
        L2_loss= 0
    loglik = torch.sum(torch.log(prob_pg))- regul*L2_loss
    loglik_org= torch.sum(torch.log(prob_pg))
    # update a
    new_a = torch.div(torch.sum(-1 * post_zg1), torch.sum(post_zg1 * torch.log(pg)))
    return loglik_org, loglik, new_a, post_zg1
    
def mstep(model, optimizer, loglike, new_a):
    optimizer.zero_grad()
    loss = -1*loglike
    loss.backward()
    optimizer.step()
    torch.clamp(model.alpha.data, min= 0.000001, max= 1, out= model.alpha.data)
    
    
def partition_data(num_genes, perc):
    train_size = int(num_genes*perc)
    all_inds = range(0, num_genes)
    train_inds = random.sample(range(0, num_genes), train_size)
    test_inds = list(set(all_inds)^set(train_inds))
    return (train_inds, test_inds)