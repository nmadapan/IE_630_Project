close all;
clear;
clc;

dlgtitle = 'IMRT Problem Variables';
prompt = {'Energy Flux (F - J/s/inch^2):', ...
          'Absorption Coeff. of normal tissue (\alpha_n - 1/inch):', ...
          'Absorption Coeff. of organ at risk (\alpha_o - 1/inch):', ...
          'Absorption Coeff. of tumor (\alpha_t - 1/inch):', ...
          'Density of normal tissue (\mu_n - g/inch^3):', ...
          'Density of organ at risk (\mu_o - g/inch^3):', ...
          'Density of of tumor (\mu_t - g/inch^3):', ...
          'Safe dosage of normal tissue (s_n - cGr):', ...
          'Safe dosage of organ at risk (s_o - cGr):', ...
          'Min. dosage of tumor (d_T - cGr):'
          };
dims = [1 50];
definput = {'10', ... % F
            '0.1', ... % \alpha_n
            '0.15', ... % \alpha_o
            '0.7', ... % \alpha_t
            '17.2', ... % \mu_n
            '20.1', ... % \mu_o
            '27', ... % \mu_t
            '0.2', ... % s_n
            '0.1', ... % s_o
            '1.05' % d_T
            };
answer = inputdlg(prompt,dlgtitle,dims,definput);
for idx = 1 : numel(answer)
    answer{idx} = str2num(answer{idx});
end
F = answer{1};
alpha_n = answer{2};
alpha_o = answer{3};
alpha_t = answer{4};
mu_n = answer{5};
mu_o = answer{6};
mu_t = answer{7};
s_n = answer{8};
s_o = answer{9};
d_T = answer{10};