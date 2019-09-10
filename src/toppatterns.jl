function toppatterns(visiblevals,visiblefeatures; patternno=10)
    freqval = freqtab(visiblevals[:,visiblefeatures])
    freqvec = reshape(freqval,prod(size(freqval)))
    topfreqs = sort(freqvec,rev=true)[1:patternno]
    toppatterns = map(index -> Float64[index.I...] .- 1,findall(freqval .>= topfreqs[end]))
end
