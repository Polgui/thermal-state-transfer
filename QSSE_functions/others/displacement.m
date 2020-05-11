function D = displacement( alpha,dimB )

    a=diag(sqrt(1:dimB-1),1);
    D=expm(alpha*a'-alpha'*a);

end

