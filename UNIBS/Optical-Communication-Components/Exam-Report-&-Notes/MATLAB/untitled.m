%%
clc
clear all

no_of_samples_1 = 1000;

V_1 = 6;
V_2 = 2;

X = linspace(-10, 10, no_of_samples_1);
X_alt = linspace(V_1, 0, no_of_samples_1);

X_disp_relation_LP01 = [(- X[0] * jv(-1, X[0])) / jv(0, X[0])];

for i=1:1:no_of_samples_1
    X_disp_relation_LP01.append( ((-X[i] * jv(-1, X[i])) / jv(0, X[i])) )
    % Y_disp_relation_LP01.append(np.sqrt(V_1**2 - X_disp_relation_LP01[i] ** 2))


Y_disp_relation_LP01 = [(X_alt[-1] * kn(-1, X_alt[-1])) / kn(0, X_alt[-1])]
%  [ np.sqrt(V_1**2 - X_disp_relation_LP01[0]**2) ]

for i in range(no_of_samples_1-1, 0, -1):
    Y_disp_relation_LP01.append( ((X_alt[i] * kn(-1, X_alt[i])) / kn(0, X_alt[i])) )


plt.figure()

plt.plot(X, X_disp_relation_LP01)
plt.plot(X_alt, Y_disp_relation_LP01, color='orange')

% label='$J_0(x)$'


plt.xlim(0, 7)
plt.ylim(-10, 10)
% plt.legend()
plt.xlabel('Don\'t Know Yet')
plt.ylabel('Don\'t Know Yet')
plt.title('Don\'t Know Yet')
% plt.savefig('../../../Graphs/2d-charge-distribution-LP01-mode.pgf', bbox_inches='tight', pad_inches=0)

plt.show()
