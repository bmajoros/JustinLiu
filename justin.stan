functions {
   // This function can be used via: p~betaModeConc(mode,concentration);
   real betaModeConc_lpdf(real parm,real m,real c) {
      return beta_lpdf(parm|m*(c-2)+1, (1-m)*(c-2)+1);
   }  
}

data {
   int<lower=0> N_RNA;      // number of RNA replicates
   int<lower=0> k[N_RNA];   // RNA alt read counts
   int<lower=0> m[N_RNA];   // RNA ref read counts
}

parameters {
   real<lower=0,upper=1> q; // overall allele freq in RNA
   real<lower=0,upper=1> qi[N_RNA]; // alt allele freqs in RNA reps
   real<lower=2> c; // concentration parameter of beta prior for qi
   real<lower=0> s; // variance parameter of lognormal prior for theta
}

model {
   c ~ gamma(1.1, 0.0005); // concentration parameter for prior on qi
   for(i in 1:N_RNA) {
      qi[i] ~ betaModeConc(q,c);
      k[i] ~ binomial(k[i]+m[i],qi[i]);
   }
}



