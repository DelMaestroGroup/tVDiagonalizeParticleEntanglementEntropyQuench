"""
Constructs a V_array with nVals elements linearly spaced between Vstart and Vend.
The array is made up of two subarrays containing positive and non-positive values
ordered with the smallest absolute values first.
"""
function lin_range(Vstart::Float64,Vend::Float64,nVals::Int64)  
    
    return collect(Vstart:abs(Vend-Vstart)/(nVals-1):Vend)

end


"""
Constructs a V_array with nVals elements logarithmically spaced between Vstart and Vend.
The array is made up of two subarrays containing positive and non-positive values
ordered with the smallest absolute values first.
"""
function log_range(Vstart::Float64,Vend::Float64,nVals::Int64) 

    if Vstart <= 0
        if Vend > 0
            v_neg = [-1.0*10^v for v in log(10,1e-1):abs(log(10,1e-1)-log(10,abs(Vstart)))/((nVals-1)/2):log(10,abs(Vstart))] 
            v_pos = [10^v for v in log(10,1e-1):abs(log(10,1e-1)-log(10,Vend))/((nVals-1)/2):log(10,Vend)] 
            
            return [reverse(v_neg);v_pos]
        else 
            return reverse([-1.0*10^v for v in  log(10,abs(Vend)):abs(log(10,abs(Vend))-log(10,abs(Vstart)))/(nVals-1):log(10,abs(Vstart))])
        end 
    end 

    return [10^v for v in log(10,Vstart):abs(log(10,Vend)-log(10,Vstart))/(nVals-1):log(10,Vend)]

end
 