//Performance should be able to be increased quadratically.
function OPFindRealEigenvector(matr,delta,limit){
    mat = cloneMat(matr)
    var vecs = []
    //
        var guess = []
               for(let m=0;m<mat.length;m++){
            guess[m]=z(Math.random()*2-1,Math.random()*2-1)
        }
        for(let n=0;n<limit;n++){
        guess = normalise(guess)
        nguess = cloneMat([guess])[0]
        for(let m=0;m<mat.length;m++){
        for(let o=0;o<mat[0].length;o++){
            var elem = cmult(guess[m],cconj(guess[o]))
            nguess[m] = csub(nguess[m],cmult(csub(mat[m][o],elem),cconj(guess[o])))
            nguess[o] = cconj(csub(cconj(nguess[o]),cmult(csub(mat[m][o],elem),guess[m])))
        }
    }
    guess=cloneMat([nguess])[0]
    }
    return normalise(guess)
}