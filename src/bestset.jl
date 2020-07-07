function expandvarset(startset,x,h=nothing,significance=1.0)
   newsize = length(startset)+1

   varno = (h == nothing) ? newsize : newsize + 1
   curx = LogLinearModels.LevelData(Matrix{Float64}(undef,size(x.data,1),varno))
   if length(startset) > 0
      curx.data[:,((h == nothing) ? 1 : 2):(end-1)] .= x.data[:,startset]
   end

   if h == nothing
      m0constr = [collect(range(1,length=varno-1)),varno]
      curx.levelno[1:(varno-1)] .= x.levelno[startset]
   else
      m0constr = [[1,collect(range(2,length=varno-2,step=1))...],varno]
      curx.levelno[2:(varno-1)] .= x.levelno[startset]
   end

   bestpartner = -1
   bestdiff = 0.0
	pval = 1.0
   for i=1:size(x.data,2)
      if i in startset
         continue
      end

      curx.data[:,varno] .= x.data[:,i]
      curx.levelno[varno] = x.levelno[i]
      curset = [startset...,i]

      for j=1:(h == nothing ? 1 : size(h.data,2))
         if h != nothing
            curx.data[:,1] .= h.data[:,j]
            curx.levelno[1] = h.levelno[j]
         end

         freqval = LogLinearModels.freqtab(curx,fillzeros=true)

         curdiff = LogLinearModels.gsquare(freqval,
                     LogLinearModels.ipf(freqval,m0constr,maxit=100))
         # println(curdiff)

         if bestpartner == -1 || curdiff > bestdiff
            bestpartner = i
            bestdiff = curdiff
         end
      end
   end
	
	if significance != 1.0
		gdiff_null = fill(0.0,100)
		for j = 1:length(gdiff_null)
			curx.data[:,varno] .= shuffle(x.data[:,bestpartner])
   		curx.levelno[varno] = x.levelno[i]
      	curset = [startset...,bestpartner]

	      for j=1:(h == nothing ? 1 : size(h.data,2))
	         if h != nothing
	            curx.data[:,1] .= h.data[:,j]
	            curx.levelno[1] = h.levelno[j]
	         end
	
	         freqval = LogLinearModels.freqtab(curx,fillzeros=true)
	
	         curdiff = LogLinearModels.gsquare(freqval,
	                     LogLinearModels.ipf(freqval,m0constr,maxit=100))
				gdiff_null[j] = curdiff
	      end
		# if (sum(gdiff_null.>=gdiff[i]) / 100) .<= significance
		pval = (sum(gdiff_null.>=bestdiff) / 100)
		end
	end
   pval, bestdiff, [startset...,bestpartner]
end

function bestpair(x,h=nothing,significance=1.0)
    p = size(x.data,2)

    bestpair = Int[]
    bestdiff = 0.0
	 pval = 1.0
    for i=1:p
        pval, curdiff, curpair = expandvarset([i],x,h,significance = significance)
        if curpair == [] || curdiff > bestdiff
            bestpair = curpair
            bestdiff = curdiff
        end
    end

    pval, bestdiff, bestpair
end

function bestset(setsize::Int,x::Array{Float64,2}; startsize=2,verbose=false,significance=1.0)
   bestset(setsize,LogLinearModels.LevelData(x),startsize=startsize,verbose=verbose,significance=significance)
end

function bestset(setsize::Int,x::Array{Float64,2},h::Array{Float64,2};
                 startsize=1,verbose=false,significance=.05)
   bestset(setsize,LogLinearModels.LevelData(x),LogLinearModels.LevelData(h),
           startsize=startsize,verbose=verbose,significance=significance)
end

function bestset(setsize::Int,x::LogLinearModels.LevelData,h=nothing;
                 startsize=1,verbose=false,significance=significance)

   gdiff = Vector{Float64}(undef,setsize - (startsize == 2 ? 1 : 0))
	pvals = Vector{Float64}(undef,setsize - (startsize == 2 ? 1 : 0))
   if h == nothing
      startsize = 2
   end

   if startsize == 2
      pvals[1], gdiff[1], curset = bestpair(x,h,significance=significance)
   else
      curset = Int[]
   end

   if verbose && length(curset) > 0
      println(curset)
   end

   for i=(startsize == 2 ? 2 : 1):(startsize == 2 ? setsize - 1 : setsize)
      pvals[i], gdiff[i], curset = expandvarset(curset,x,h,significance)
		
      if verbose
         println(curset)
      end
   end
   println("d")
   return gdiff,curset
end
