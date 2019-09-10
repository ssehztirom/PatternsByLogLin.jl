function kldivergence(P,Q)
    sum(broadcast((pval, qval) -> pval != 0.0 ? pval*log(pval/qval) : 0.0,P,Q))
end

function jsdivergence(P,Q)
    M = (P .+ Q) ./ 2
    0.5*kldivergence(P,M) + 0.5*kldivergence(Q,M)
end

function patterncoords(visiblevals,hiddenvals,visiblefeatures,patterns)
    topindices = map(pattern -> [all(visiblevals[i,visiblefeatures] .== pattern)
                                for i=1:size(visiblevals,1)],
                     patterns)
    topfreqs = map(indices -> sum(indices),topindices)
    tophidden = map(indices -> hiddenvals[indices,:],topindices)

    tophrelfreqval = map(tophidden) do hidden
        hfreqvals = freqtab(hidden)
        hfreqvals ./ sum(hfreqvals)
    end

    distmat = fill(0.0,length(tophidden),length(tophidden))
    for i=1:(length(tophidden)-1), j=(i+1):length(tophidden)
        distmat[i,j] = distmat[j,i] = sqrt(jsdivergence(tophrelfreqval[i],tophrelfreqval[j]))
    end

    mdscoord = MultivariateStats.classical_mds(distmat,2)

    distmat, mdscoord
end
