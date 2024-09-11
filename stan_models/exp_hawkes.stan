functions {
    matrix zero_above_diagonal(matrix input_mat){
        int N = rows(input_mat);
        matrix[N, N] output_mat;

        for (i in 1:N) {
            for (j in 1:N) {
                if (j >= i) {
                   output_mat[i, j] = 0;
                } else {
                   output_mat[i, j] = input_mat[i, j];
                }
            }
        }

        return output_mat;
    }

    vector rowsum(matrix input_mat){
        int M = cols(input_mat);
        vector[M] v = rep_vector(1.0, M);

        return input_mat * v;
    }

    vector non_negative_mask(vector input_vec){
        int N = num_elements(input_vec);
        vector[N] output_vec;

        for (i in 1:N) {
            output_vec[i] = input_vec[i] > 0;
        }

        return output_vec;
    }

    matrix self_differences(vector input_vec){
        int N = num_elements(input_vec);

        row_vector[N] rv = input_vec';

        matrix[N, N] output_mat = rep_matrix(input_vec, N) - rep_matrix(rv, N);

        return output_mat;
    }

    real log_likelihood(real mu, real alpha, real delta, vector events_list, int N, real max_T){

        // first term

        vector[N] differences_from_max_T = max_T - events_list;
        vector[N] summands = exp(-delta * differences_from_max_T) - 1;

        vector[N] within_max_T_Mask = non_negative_mask(differences_from_max_T);
        summands = summands .* within_max_T_Mask;

        real first = mu * max_T - (alpha / delta) * sum(summands);

        // second term

        matrix[N, N] differences_mat = self_differences(events_list);
        matrix[N, N] inner_sum_mat = exp(-delta * differences_mat);
        inner_sum_mat = zero_above_diagonal(inner_sum_mat);

        vector[N] term_inside_log = mu + alpha * rowsum(inner_sum_mat);

        vector[N] second_sum_terms = log(term_inside_log);

        real second = sum(second_sum_terms);

        return -first + second;
    }
}

data {
    int<lower=0> N;
    vector[N] events_list;
    real max_T;
}

parameters {
    real mu;
    real <lower=0> alpha;
    real <lower=0> delta;
}

model {
    mu ~ normal( 0.4 , 1 );
    alpha ~ normal( 0.4 , 0.1 );
    delta ~ normal( 0.5 , 1 );

    target += log_likelihood(mu, alpha, delta, events_list, N, max_T);
}

generated quantities {
  real log_lik = log_likelihood(mu, alpha, delta, events_list, N, max_T);
}
