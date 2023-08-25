
import numpy as np
from tqdm import tqdm

# A可能不是方阵，A是numpy类型的数组
# 这个函数已经测试完毕
def Jaccard_similarity(A, A1):
    # B是返回的杰卡德相似性矩阵
    B=np.zeros((A.shape[0],A.shape[0]))
    for i in tqdm(range(B.shape[0])):
        for j in range(i+1,B.shape[1]):#只做上三角部分

            if np.sum(A[i])==0 and np.sum(A1[j])==0:#如果两个药物不和任何靶点相互作用，则两个药物的相似度为0
                B[i][j]=0
            else:
                jiaoji=0
                bingji=0
                # 计算A[i]和A[j]的交集和并集
                for k in range(A.shape[1]):
                    for k1 in range(A1.shape[1]):
                        if A[i][k]==1 and A1[j][k1]==1:
                           jiaoji+=1
                           bingji+=1
                        elif A[i][k]==1 or A1[j][k1]==1:
                           bingji+=1
                B[i][j]=jiaoji/bingji
    # 此时B只是上三角矩阵，将上三角矩阵转为对称阵
    # 将主对角元素置为1
    # 因为有些药物不和任何靶点相互作用，但是自己和自己的相似度肯定是1
    row,col=np.diag_indices_from(B)
    B[row,col]=1
    B += B.T - np.diag(B.diagonal())

    # print(B.T==B)


    return B

if __name__=="__main__":
   # A=np.array([[0,0,0,1],
    #         [1,1,0,0],
     #           [0,1,0,0],
      #           [1,1,1,1],
       #         [1,0,1,0],
       #          [0,0,0,0]])

  #  A1 = np.array([[0, 1, 1, 0],
    #            [0, 0, 0, 0],
     #             [0, 0, 0, 0],
     #              [1, 1, 0, 0],
      #             [1, 1, 1, 1],
        #           [0, 0, 1, 1]])
   # A=np.array(A)
   # A1 = np.array(A1)
     # ans=Jaccard_similarity(A)
   # drug_target_interaction= Jaccard_similarity(A, A1)

# 药物-药物相互作用矩阵





    A = np.loadtxt('dataset/mat_drug_disease.txt')



    A1 = np.loadtxt('dataset/mat_protein_disease.txt')


    drug_t_interaction = np.loadtxt('dataset/mat_drug_protein.txt')

    drug_target_interaction = Jaccard_similarity(A, A1)


    #np.save('dataset/multi_similarity/drug_target_interaction.npy', drug_target_interaction)


    x1 = np.maximum(drug_target_interaction, drug_t_interaction)


    np.savetxt('dataset/multi_similarity/DTI_708_1512_MAX_DISCRETIZE.txt', x1)





#print('end')

print('B')