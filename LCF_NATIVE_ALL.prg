' LCF_NATIVE_ALL.prg -- all statistical computation is native EViews.
' Self-contained: the raw annual data and required reference constants are embedded.
' Save this PRG in a writable folder. Open in EViews and Run.
' Outputs, including LCF_NATIVE_RESULTS.wf1, are written beside this PRG.
' Full EViews runtime execution has NOT been tested in the author's Linux environment.
' Check LCF_STATUS.csv for SUCCESS before interpreting outputs as a complete run.
' Seeds use the EViews generator, so bootstrap draws differ from other software.
' Existing same-named LCF_* outputs in this directory will be replaced on rerun.
' ---------------------------------------------------------------------------
' Build: 2026-09-26-r6; progress: LCF_PROGRESS.csv; partial results: LCF_models_PARTIAL.csv.
subroutine local olsfit(matrix x, vector y, vector b, vector e, matrix gi, scalar ss, scalar rk)
  !of_n=@rows(x)
  !of_k=@columns(x)
  vector(!of_k) sc=1
  for !of_j=1 to !of_k
    sc(!of_j)=@sqrt(@sum(@epow(@columnextract(x,!of_j),2)))
    if sc(!of_j)=0 then
      sc(!of_j)=1
    endif
  next
  matrix xs=@ediv(x,@ones(!of_n,1)*@transpose(sc))
  rk=@rank(xs)
  ' Direct thin SVD: never form X'X, or call pinverse on a column matrix.
  vector singular
  matrix rightvec
  matrix leftvec=@svd(xs,singular,rightvec)
  !of_r=@rows(singular)
  if @rows(leftvec)<>!of_n or @columns(leftvec)<>!of_r or @rows(rightvec)<>!of_k or @columns(rightvec)<>!of_r then
    seterr "Unexpected EViews SVD dimensions. Send the error and LCF_PROGRESS.csv."
  endif
  matrix inverse_s=@zeros(!of_r,!of_r)
  !of_largest=!of_n
  if !of_k>!of_largest then
    !of_largest=!of_k
  endif
  !of_tol=!of_largest*@max(singular)*2.220446049250313e-16
  for !of_j=1 to !of_r
    if singular(!of_j)>!of_tol then
      inverse_s(!of_j,!of_j)=1/singular(!of_j)
    endif
  next
  matrix pi=rightvec*inverse_s*@transpose(leftvec)
  pi=@ediv(pi,sc*@ones(1,!of_n))
  b=pi*y
  gi=pi*@transpose(pi)
  e=y-x*b
  ss=@sum(@epow(e,2))
endsub

subroutine local design(vector y, matrix x, vector ord, scalar det, scalar hold, scalar ecm, matrix w, vector v, scalar nd)
  !ds_T=@rows(y)
  !ds_n=!ds_T-hold
  nd=1+(det=2)+(det=4)+(det=3)+(det=4)
  !ds_p=ord(1)
  !ds_k=nd+!ds_p+5
  for !ds_j=2 to 6
    !ds_k=!ds_k+ord(!ds_j)
  next
  w=@ones(!ds_n,!ds_k)
  vector piece
  !ds_col=1
  if det=2 or det=4 then
    !ds_col=!ds_col+1
    for !ds_i=1 to !ds_n
      w(!ds_i,!ds_col)=(hold+!ds_i-1)/60
    next
  endif
  if det=3 or det=4 then
    !ds_col=!ds_col+1
    for !ds_i=1 to !ds_n
      w(!ds_i,!ds_col)=(hold+!ds_i+1964>=1979)
    next
  endif
  v=@subextract(y,hold+1,1,!ds_T,1)
  if ecm=0 then
    for !ds_l=1 to !ds_p
      !ds_col=!ds_col+1
      piece=@subextract(y,hold+1-!ds_l,1,!ds_T-!ds_l,1)
      colplace(w,piece,!ds_col)
    next
    for !ds_j=1 to 5
      for !ds_l=0 to ord(!ds_j+1)
        !ds_col=!ds_col+1
        piece=@subextract(x,hold+1-!ds_l,!ds_j,!ds_T-!ds_l,!ds_j)
        colplace(w,piece,!ds_col)
      next
    next
  else
    piece=@subextract(y,hold,1,!ds_T-1,1)
    v=v-piece
    !ds_col=!ds_col+1
    colplace(w,piece,!ds_col)
    for !ds_j=1 to 5
      !ds_col=!ds_col+1
      !ds_l=1
      if ord(!ds_j+1)=0 then
        ' Exact q=0 tied level/difference coefficient: use current x.
        !ds_l=0
      endif
      piece=@subextract(x,hold+1-!ds_l,!ds_j,!ds_T-!ds_l,!ds_j)
      colplace(w,piece,!ds_col)
    next
    if !ds_p>1 then
      for !ds_l=1 to !ds_p-1
        !ds_col=!ds_col+1
        piece=@subextract(y,hold+1-!ds_l,1,!ds_T-!ds_l,1)-@subextract(y,hold-!ds_l,1,!ds_T-!ds_l-1,1)
        colplace(w,piece,!ds_col)
      next
    endif
    for !ds_j=1 to 5
      if ord(!ds_j+1)>0 then
        for !ds_l=0 to ord(!ds_j+1)-1
          !ds_col=!ds_col+1
          piece=@subextract(x,hold+1-!ds_l,!ds_j,!ds_T-!ds_l,!ds_j)-@subextract(x,hold-!ds_l,!ds_j,!ds_T-!ds_l-1,!ds_j)
          colplace(w,piece,!ds_col)
        next
      endif
    next
  endif
endsub

subroutine local export_table(matrix a, string %headers, string %filename, table tab)
  !et_n=@rows(a)
  !et_k=@columns(a)
  for !et_j=1 to !et_k
    tab(1,!et_j)=@word(%headers,!et_j)
  next
  for !et_i=1 to !et_n
    for !et_j=1 to !et_k
      tab(!et_i+1,!et_j)=a(!et_i,!et_j)
    next
  next
  tab.setformat(@all) g.12
  !et_lastrow=!et_n+1
  tab.save(t=csv,f,r=R1C1:R{!et_lastrow}C{!et_k}) %filename
endsub

subroutine local lagsearch(vector y, matrix x, scalar cap, scalar det, scalar hold, scalar qmin, vector obic, vector oaicc, matrix candidates)
  !ls_count=cap*(cap-qmin+1)^5
  candidates=@zeros(!ls_count,11)
  vector(6) ord
  matrix w
  vector v
  vector b
  vector e
  matrix gi
  scalar ss
  scalar rk
  scalar nd
  !ls_bestb=1e100
  !ls_besta=1e100
  !ls_row=0
  !ls_done=0
  for !ls_p=1 to cap
    for !ls_q1=qmin to cap
      for !ls_q2=qmin to cap
        for !ls_q3=qmin to cap
          for !ls_q4=qmin to cap
            for !ls_q5=qmin to cap
              if !ls_done=0 or !ls_done=250*@floor(!ls_done/250) then
                %ls_msg="Lag search L="+@str(cap)+" det="+@str(det)+" qmin="+@str(qmin)+" hold="+@str(hold)
                call progress(%ls_msg,!ls_done,!ls_count,0)
              endif
              !ls_done=!ls_done+1
              ord(1)=!ls_p
              ord(2)=!ls_q1
              ord(3)=!ls_q2
              ord(4)=!ls_q3
              ord(5)=!ls_q4
              ord(6)=!ls_q5
              call design(y,x,ord,det,hold,0,w,v,nd)
              call olsfit(w,v,b,e,gi,ss,rk)
              !ls_n=@rows(w)
              !ls_k=@columns(w)
              if rk=!ls_k and !ls_n-!ls_k>=12 then
                !ls_row=!ls_row+1
                !ls_kv=!ls_k+1
                !ls_ll=-!ls_n/2*(log(2*3.141592653589793)+1+log(ss/!ls_n))
                !ls_bic=-2*!ls_ll+!ls_kv*log(!ls_n)
                !ls_aicc=-2*!ls_ll+2*!ls_kv+2*!ls_kv*(!ls_kv+1)/(!ls_n-!ls_kv-1)
                for !ls_j=1 to 6
                  candidates(!ls_row,!ls_j)=ord(!ls_j)
                next
                candidates(!ls_row,7)=!ls_n
                candidates(!ls_row,8)=!ls_k
                candidates(!ls_row,9)=!ls_ll
                candidates(!ls_row,10)=!ls_bic
                candidates(!ls_row,11)=!ls_aicc
                if !ls_bic<!ls_bestb then
                  !ls_bestb=!ls_bic
                  obic=ord
                endif
                if !ls_aicc<!ls_besta then
                  !ls_besta=!ls_aicc
                  oaicc=ord
                endif
              endif
            next
          next
        next
      next
    next
  next
  if !ls_row=0 then
    seterr "No admissible full-rank lag candidate."
  endif
  call progress("Lag search completed",!ls_done,!ls_count,0)
  candidates=@subextract(candidates,1,1,!ls_row,11)
endsub

subroutine local teststats(vector b, matrix gi, scalar ss, scalar n, scalar nd, vector st)
  !ts_k=@rows(b)
  matrix cv=gi*(ss/(n-!ts_k))
  vector bl=@subextract(b,nd+1,1,nd+6,1)
  matrix vl=@subextract(cv,nd+1,nd+1,nd+6,nd+6)
  matrix z=@transpose(bl)*@inverse(vl)*bl
  st=@zeros(3)
  st(1)=z(1,1)/6
  st(2)=b(nd+1)/@sqrt(cv(nd+1,nd+1))
  bl=@subextract(b,nd+2,1,nd+6,1)
  vl=@subextract(cv,nd+2,nd+2,nd+6,nd+6)
  z=@transpose(bl)*@inverse(vl)*bl
  st(3)=z(1,1)/5
endsub

subroutine local quantile7(vector draws, scalar prob, scalar value)
  vector sorted=@sort(draws)
  !qt_pos=1+(@rows(draws)-1)*prob
  !qt_lo=@floor(!qt_pos)
  !qt_hi=!qt_lo+1
  if !qt_hi>@rows(draws) then
    !qt_hi=@rows(draws)
  endif
  value=sorted(!qt_lo)+(!qt_pos-!qt_lo)*(sorted(!qt_hi)-sorted(!qt_lo))
endsub

subroutine local recurse_y(vector y, matrix x, vector ord, scalar det, scalar hold, matrix w, scalar nd, vector b, vector eps, vector ys)
  ys=y
  vector fixed_b=b
  fixed_b(nd+1)=0
  if ord(1)>1 then
    for !ry_l=1 to ord(1)-1
      fixed_b(nd+6+!ry_l)=0
    next
  endif
  vector fixed_part=w*fixed_b
  for !ry_i=1 to @rows(w)
    !ry_t=hold+!ry_i
    !ry_mu=fixed_part(!ry_i)+b(nd+1)*ys(!ry_t-1)
    if ord(1)>1 then
      for !ry_l=1 to ord(1)-1
        !ry_mu=!ry_mu+b(nd+6+!ry_l)*(ys(!ry_t-!ry_l)-ys(!ry_t-!ry_l-1))
      next
    endif
    ys(!ry_t)=ys(!ry_t-1)+!ry_mu+eps(!ry_i)
  next
endsub

subroutine local bounds_boot(vector y, matrix x, vector ord, scalar det, scalar hold, scalar reps, scalar seed, scalar iid, scalar model_id, matrix answer)
  !bb_rng_seed=seed
  rndseed {!bb_rng_seed}
  matrix w
  vector v
  scalar nd
  call design(y,x,ord,det,hold,1,w,v,nd)
  !bb_n=@rows(w)
  !bb_k=@columns(w)
  vector b
  vector e
  matrix gi
  scalar ss
  scalar rk
  vector st
  call olsfit(w,v,b,e,gi,ss,rk)
  call teststats(b,gi,ss,!bb_n,nd,st)
  answer=@zeros(3,4)
  vector(!bb_k) keep
  vector br
  vector er
  matrix gr
  scalar sr
  scalar rr
  matrix wr
  vector(!bb_n) innov
  vector(!bb_n) eps
  vector(!bb_k) bfull
  vector ys
  matrix ws
  vector vs
  vector bs
  vector es
  matrix gs
  scalar sss
  scalar rks
  scalar nds
  vector sts
  vector(reps) draws
  scalar crit
  for !bb_test=1 to 3
    !bb_nr=0
    for !bb_j=1 to !bb_k
      !bb_zero=0
      if !bb_test=1 and !bb_j>=nd+1 and !bb_j<=nd+6 then
        !bb_zero=1
      endif
      if !bb_test=2 and !bb_j=nd+1 then
        !bb_zero=1
      endif
      if !bb_test=3 and !bb_j>=nd+2 and !bb_j<=nd+6 then
        !bb_zero=1
      endif
      if !bb_zero=0 then
        !bb_nr=!bb_nr+1
        keep(!bb_nr)=!bb_j
      endif
    next
    wr=@zeros(!bb_n,!bb_nr)
    for !bb_j=1 to !bb_nr
      for !bb_i=1 to !bb_n
        wr(!bb_i,!bb_j)=w(!bb_i,keep(!bb_j))
      next
    next
    call olsfit(wr,v,br,er,gr,sr,rr)
    bfull=0
    for !bb_j=1 to !bb_nr
      bfull(keep(!bb_j))=br(!bb_j)
    next
    for !bb_i=1 to !bb_n
      !bb_hat=0
      for !bb_j=1 to !bb_nr
        for !bb_l=1 to !bb_nr
          !bb_hat=!bb_hat+wr(!bb_i,!bb_j)*gr(!bb_j,!bb_l)*wr(!bb_i,!bb_l)
        next
      next
      !bb_d=1-!bb_hat
      if !bb_d<0.05 then
        !bb_d=0.05
      endif
      innov(!bb_i)=er(!bb_i)/@sqrt(!bb_d)
      if iid=1 then
        innov(!bb_i)=(er(!bb_i)-@mean(er))*@sqrt(!bb_n/(!bb_n-!bb_nr))
      endif
    next
    !bb_tail=0
    for !bb_rep=1 to reps
      if !bb_rep=1 or (!bb_rep-1)=50*@floor((!bb_rep-1)/50) then
        %bb_msg="Bounds test="+@str(!bb_test)+" (1=F 2=t 3=Fx); iid="+@str(iid)
        call progress(%bb_msg,!bb_rep-1,reps,model_id)
      endif
      for !bb_i=1 to !bb_n
        if iid=0 then
          eps(!bb_i)=innov(!bb_i)*(2*(@runif(0,1)>=0.5)-1)
        else
          !bb_pick=1+@floor(@runif(0,1)*!bb_n)
          eps(!bb_i)=innov(!bb_pick)
        endif
      next
      call recurse_y(y,x,ord,det,hold,w,nd,bfull,eps,ys)
      call refresh_ecm(ys,ord,hold,w,nd,ws,vs)
      nds=nd
      call olsfit(ws,vs,bs,es,gs,sss,rks)
      if rks<>@columns(ws) then
        seterr "Rank failure in bootstrap. No invalid draw was silently discarded."
      endif
      call teststats(bs,gs,sss,!bb_n,nds,sts)
      draws(!bb_rep)=sts(!bb_test)
      if !bb_test=2 then
        !bb_tail=!bb_tail+(draws(!bb_rep)<=st(!bb_test))
      else
        !bb_tail=!bb_tail+(draws(!bb_rep)>=st(!bb_test))
      endif
    next
    !bb_prob=0.95
    if !bb_test=2 then
      !bb_prob=0.05
    endif
    call quantile7(draws,!bb_prob,crit)
    answer(!bb_test,1)=st(!bb_test)
    answer(!bb_test,2)=(1+!bb_tail)/(reps+1)
    answer(!bb_test,3)=crit
    answer(!bb_test,4)=@sqrt(answer(!bb_test,2)*(1-answer(!bb_test,2))/(reps+1))
  next
endsub

subroutine local covariances(matrix w, vector e, matrix gi, scalar ss, matrix olscov, matrix hac, vector lev, vector cook)
  !cv_n=@rows(w)
  !cv_k=@columns(w)
  olscov=gi*(ss/(!cv_n-!cv_k))
  matrix scores=w
  lev=@zeros(!cv_n)
  cook=@zeros(!cv_n)
  for !cv_i=1 to !cv_n
    for !cv_j=1 to !cv_k
      scores(!cv_i,!cv_j)=w(!cv_i,!cv_j)*e(!cv_i)
      for !cv_l=1 to !cv_k
        lev(!cv_i)=lev(!cv_i)+w(!cv_i,!cv_j)*gi(!cv_j,!cv_l)*w(!cv_i,!cv_l)
      next
    next
    cook(!cv_i)=e(!cv_i)^2/(ss/(!cv_n-!cv_k))/!cv_k*lev(!cv_i)/(1-lev(!cv_i))^2
  next
  matrix meat=@transpose(scores)*scores
  matrix c1
  for !cv_lag=1 to 2
    c1=@transpose(@subextract(scores,!cv_lag+1,1,!cv_n,!cv_k))*@subextract(scores,1,1,!cv_n-!cv_lag,!cv_k)
    meat=meat+(1-!cv_lag/3)*(c1+@transpose(c1))
  next
  hac=gi*meat*gi*!cv_n/(!cv_n-!cv_k)
endsub

subroutine local long_run(vector b, matrix cv, scalar nd, scalar df, matrix lr)
  !lr_k=@rows(b)
  lr=@zeros(8,5)
  vector(!lr_k) grad
  matrix vv
  !lr_crit=@qtdist(0.975,df)
  for !lr_j=1 to 8
    grad=0
    if !lr_j<=5 then
      !lr_i=nd+1+!lr_j
      !lr_num=b(!lr_i)
      grad(!lr_i)=-1/b(nd+1)
    else
      !lr_s=!lr_j-7
      !lr_num=b(nd+4)+!lr_s*b(nd+6)
      grad(nd+4)=-1/b(nd+1)
      grad(nd+6)=-!lr_s/b(nd+1)
    endif
    grad(nd+1)=!lr_num/b(nd+1)^2
    vv=@transpose(grad)*cv*grad
    lr(!lr_j,1)=-!lr_num/b(nd+1)
    lr(!lr_j,2)=@sqrt(vv(1,1))
    lr(!lr_j,3)=lr(!lr_j,1)-!lr_crit*lr(!lr_j,2)
    lr(!lr_j,4)=lr(!lr_j,1)+!lr_crit*lr(!lr_j,2)
    lr(!lr_j,5)=2*(1-@ctdist(@abs(lr(!lr_j,1)/lr(!lr_j,2)),df))
  next
endsub

subroutine local diagnostics(matrix w, vector v, vector b, vector e, matrix gi, scalar ss, matrix a, vector ya, vector ba, vector diagout)
  !dg_n=@rows(w)
  !dg_k=@columns(w)
  diagout=@zeros(24)
  matrix aux
  vector target
  vector bb
  vector ee
  matrix gg
  scalar sse
  scalar ranka
  !dg_sst=@sum(@epow(e-@mean(e),2))
  !dg_lb=0
  !dg_slot=0
  vector(3) lags
  lags(1)=1
  lags(2)=2
  lags(3)=4
  for !dg_h=1 to 3
    !dg_l=lags(!dg_h)
    aux=@zeros(!dg_n,!dg_k+!dg_l)
    for !dg_i=1 to !dg_n
      for !dg_j=1 to !dg_k
        aux(!dg_i,!dg_j)=w(!dg_i,!dg_j)
      next
      for !dg_j=1 to !dg_l
        if !dg_i>!dg_j then
          aux(!dg_i,!dg_k+!dg_j)=e(!dg_i-!dg_j)
        endif
      next
    next
    call olsfit(aux,e,bb,ee,gg,sse,ranka)
    diagout(!dg_h)=1-@cchisq(!dg_n*(1-sse/!dg_sst),!dg_l)
    !dg_f=((ss-sse)/!dg_l)/(sse/(!dg_n-!dg_k-!dg_l))
    diagout(3+!dg_h)=1-@cfdist(!dg_f,!dg_l,!dg_n-!dg_k-!dg_l)
    !dg_lb=0
    for !dg_lag=1 to !dg_l
      !dg_ac=0
      for !dg_i=!dg_lag+1 to !dg_n
        !dg_ac=!dg_ac+(e(!dg_i)-@mean(e))*(e(!dg_i-!dg_lag)-@mean(e))
      next
      !dg_ac=!dg_ac/!dg_sst
      !dg_lb=!dg_lb+!dg_ac^2/(!dg_n-!dg_lag)
    next
    diagout(6+!dg_h)=1-@cchisq(!dg_n*(!dg_n+2)*!dg_lb,!dg_l)
    aux=@ones(!dg_n-!dg_l,!dg_l+1)
    target=@zeros(!dg_n-!dg_l)
    for !dg_i=1 to !dg_n-!dg_l
      target(!dg_i)=e(!dg_i+!dg_l)^2
      for !dg_j=1 to !dg_l
        aux(!dg_i,!dg_j+1)=e(!dg_i+!dg_l-!dg_j)^2
      next
    next
    call olsfit(aux,target,bb,ee,gg,sse,ranka)
    !dg_r2=1-sse/@sum(@epow(target-@mean(target),2))
    diagout(9+!dg_h)=1-@cchisq((!dg_n-!dg_l)*!dg_r2,!dg_l)
  next
  target=@epow(e,2)
  call olsfit(w,target,bb,ee,gg,sse,ranka)
  !dg_r2=1-sse/@sum(@epow(target-@mean(target),2))
  diagout(13)=1-@cchisq(!dg_n*!dg_r2,!dg_k-1)
  !dg_m2=@mean(@epow(e-@mean(e),2))
  !dg_sk=@mean(@epow(e-@mean(e),3))/!dg_m2^1.5
  !dg_ku=@mean(@epow(e-@mean(e),4))/!dg_m2^2
  diagout(14)=1-@cchisq(!dg_n/6*(!dg_sk^2+(!dg_ku-3)^2/4),2)
  !dg_dw=0
  for !dg_i=2 to !dg_n
    !dg_dw=!dg_dw+(e(!dg_i)-e(!dg_i-1))^2
  next
  diagout(15)=!dg_dw/ss
  vector fitted=w*b
  aux=@hcat(w,@epow(fitted,2))
  call olsfit(aux,v,bb,ee,gg,sse,ranka)
  diagout(16)=1-@cfdist(((ss-sse)/(sse/(!dg_n-!dg_k-1))),1,!dg_n-!dg_k-1)
  fitted=a*ba
  aux=@hcat(a,@epow(fitted,2))
  call olsfit(aux,ya,bb,ee,gg,sse,ranka)
  diagout(17)=1-@cfdist(((ss-sse)/(sse/(!dg_n-!dg_k-1))),1,!dg_n-!dg_k-1)
  !dg_sum=0
  !dg_sup=0
  for !dg_i=1 to !dg_n
    !dg_sum=!dg_sum+e(!dg_i)
    if @abs(!dg_sum)>!dg_sup then
      !dg_sup=@abs(!dg_sum)
    endif
  next
  !dg_stat=!dg_sup/@sqrt(ss*!dg_n/(!dg_n-!dg_k))
  !dg_prob=0
  for !dg_j=1 to 1000
    !dg_term=2*(-1)^(!dg_j-1)*exp(-2*!dg_j^2*!dg_stat^2)
    !dg_prob=!dg_prob+!dg_term
    if @abs(!dg_term)<1e-15 then
      exitloop
    endif
  next
  diagout(18)=!dg_stat
  diagout(19)=!dg_prob
  aux=@zeros(!dg_n,!dg_k*(!dg_k+1)/2)
  !dg_col=0
  for !dg_j=1 to !dg_k
    for !dg_l=!dg_j to !dg_k
      !dg_col=!dg_col+1
      for !dg_i=1 to !dg_n
        aux(!dg_i,!dg_col)=w(!dg_i,!dg_j)*w(!dg_i,!dg_l)
      next
    next
  next
  !dg_wr=@rank(aux)
  diagout(20)=!dg_wr
  diagout(21)=!dg_n-!dg_wr
  diagout(22)=NA
  if !dg_n-!dg_wr>0 then
    target=@epow(e,2)
    call olsfit(aux,target,bb,ee,gg,sse,ranka)
    !dg_r2=1-sse/@sum(@epow(target-@mean(target),2))
    diagout(22)=1-@cchisq(!dg_n*!dg_r2,!dg_wr-1)
  endif
  vector sv
  matrix vt
  matrix uu=@svd(w,sv,vt)
  diagout(23)=@max(sv)/@min(sv)
  diagout(24)=!dg_n-!dg_k
endsub

subroutine local ar_modulus(vector phi, scalar modulus)
  ' Durand-Kerner iteration on lambda^p-phi1*lambda^(p-1)-...-phip.
  !ar_p=@rows(phi)
  if !ar_p=1 then
    modulus=@abs(phi(1))
    return
  endif
  vector(!ar_p) re
  vector(!ar_p) im
  !ar_radius=1+@max(@abs(phi))
  for !ar_j=1 to !ar_p
    !ar_angle=2*3.141592653589793*(!ar_j-1)/!ar_p+0.19
    re(!ar_j)=!ar_radius*cos(!ar_angle)
    im(!ar_j)=!ar_radius*sin(!ar_angle)
  next
  for !ar_it=1 to 1000
    !ar_move=0
    for !ar_j=1 to !ar_p
      !ar_pr=1
      !ar_pi=0
      !ar_dr=1
      !ar_di=0
      for !ar_l=1 to !ar_p
        !ar_new=!ar_pr*re(!ar_j)-!ar_pi*im(!ar_j)-phi(!ar_l)
        !ar_pi=!ar_pr*im(!ar_j)+!ar_pi*re(!ar_j)
        !ar_pr=!ar_new
        if !ar_l<>!ar_j then
          !ar_rr=re(!ar_j)-re(!ar_l)
          !ar_ii=im(!ar_j)-im(!ar_l)
          !ar_new=!ar_dr*!ar_rr-!ar_di*!ar_ii
          !ar_di=!ar_dr*!ar_ii+!ar_di*!ar_rr
          !ar_dr=!ar_new
        endif
      next
      !ar_denom=!ar_dr^2+!ar_di^2
      if !ar_denom<1e-30 then
        seterr "AR-root iteration became singular."
      endif
      !ar_delta_r=(!ar_pr*!ar_dr+!ar_pi*!ar_di)/!ar_denom
      !ar_delta_i=(!ar_pi*!ar_dr-!ar_pr*!ar_di)/!ar_denom
      re(!ar_j)=re(!ar_j)-!ar_delta_r
      im(!ar_j)=im(!ar_j)-!ar_delta_i
      !ar_move=!ar_move+!ar_delta_r^2+!ar_delta_i^2
    next
    if !ar_move<1e-24 then
      exitloop
    endif
  next
  if !ar_move>=1e-24 then
    seterr "AR-root iteration did not converge."
  endif
  modulus=0
  for !ar_j=1 to !ar_p
    !ar_abs=@sqrt(re(!ar_j)^2+im(!ar_j)^2)
    if !ar_abs>modulus then
      modulus=!ar_abs
    endif
  next
endsub

subroutine local supf(matrix w, vector v, scalar value, scalar cut)
  !sf_n=@rows(w)
  !sf_k=@columns(w)
  vector b
  vector e
  matrix gi
  scalar ss
  scalar rk
  call olsfit(w,v,b,e,gi,ss,rk)
  !sf_sse=ss
  !sf_lo=@floor(0.2*!sf_n)
  if !sf_lo<!sf_k+5 then
    !sf_lo=!sf_k+5
  endif
  !sf_hi=@floor(0.8*!sf_n)
  if !sf_hi>!sf_n-!sf_k-4 then
    !sf_hi=!sf_n-!sf_k-4
  endif
  matrix split
  matrix w2
  value=NA
  cut=-1
  if !sf_lo>=!sf_hi then
    return
  endif
  for !sf_c=!sf_lo to !sf_hi-1
    w2=w
    matrix zeroblock=@zeros(!sf_c,!sf_k)
    matplace(w2,zeroblock,1,1)
    split=@hcat(w,w2)
    if @rank(split)=2*!sf_k then
      call olsfit(split,v,b,e,gi,ss,rk)
      !sf_val=((!sf_sse-ss)/!sf_k)/(ss/(!sf_n-2*!sf_k))
      if cut=-1 then
        value=!sf_val
        cut=!sf_c
      else
        if !sf_val>value then
          value=!sf_val
          cut=!sf_c
        endif
      endif
    endif
  next
endsub

subroutine local recursive_paths(matrix w, vector v, scalar first, scalar csmax, scalar sqmax, matrix path)
  !rp_n=@rows(w)
  !rp_k=@columns(w)
  !rp_m=!rp_n-first
  vector(!rp_m) rec
  matrix wa
  vector va
  vector bb
  vector ee
  matrix gg
  scalar ss
  scalar rk
  for !rp_a=first to !rp_n-1
    wa=@subextract(w,1,1,!rp_a,!rp_k)
    va=@subextract(v,1,1,!rp_a,1)
    call olsfit(wa,va,bb,ee,gg,ss,rk)
    matrix nextrow=@subextract(w,!rp_a+1,1,!rp_a+1,!rp_k)
    matrix prediction=nextrow*bb
    matrix predvariance=nextrow*gg*@transpose(nextrow)
    !rp_fit=prediction(1,1)
    !rp_var=1+predvariance(1,1)
    rec(!rp_a-first+1)=(v(!rp_a+1)-!rp_fit)/@sqrt(!rp_var)
  next
  !rp_sd=@sqrt(@sum(@epow(rec-@mean(rec),2))/(!rp_m-1))
  !rp_ss=@sum(@epow(rec,2))
  path=@zeros(!rp_m,2)
  !rp_cum=0
  !rp_sq=0
  csmax=0
  sqmax=0
  for !rp_i=1 to !rp_m
    !rp_cum=!rp_cum+rec(!rp_i)
    !rp_sq=!rp_sq+rec(!rp_i)^2
    path(!rp_i,1)=!rp_cum/(!rp_sd*@sqrt(!rp_m))
    path(!rp_i,2)=!rp_sq/!rp_ss-!rp_i/!rp_m
    if @abs(path(!rp_i,1))>csmax then
      csmax=@abs(path(!rp_i,1))
    endif
    if @abs(path(!rp_i,2))>sqmax then
      sqmax=@abs(path(!rp_i,2))
    endif
  next
endsub

subroutine local stability_boot(vector y, matrix x, vector ord, scalar det, scalar hold, scalar reps, scalar model_id, vector answer, matrix observedpath)
  rndseed 64217
  matrix w
  vector v
  scalar nd
  call design(y,x,ord,det,hold,1,w,v,nd)
  !sb_n=@rows(w)
  !sb_k=@columns(w)
  scalar first=-1
  for !sb_a=!sb_k+2 to !sb_n-6
    if @rank(@subextract(w,1,1,!sb_a,!sb_k))=!sb_k then
      first=!sb_a
      exitloop
    endif
  next
  if first<0 then
    seterr "No admissible full-rank initial window for stability."
  endif
  vector b
  vector e
  matrix gi
  scalar ss
  scalar rk
  call olsfit(w,v,b,e,gi,ss,rk)
  scalar fs
  scalar cut
  scalar cm
  scalar qm
  call supf(w,v,fs,cut)
  call recursive_paths(w,v,first,cm,qm,observedpath)
  answer=@zeros(10)
  answer(1)=fs
  answer(2)=1965+hold+cut
  if cut<0 then
    answer(2)=NA
  endif
  answer(4)=cm
  answer(6)=qm
  answer(8)=1965+hold+first
  vector(!sb_n) innov
  vector(!sb_n) eps
  for !sb_i=1 to !sb_n
    !sb_hat=0
    for !sb_j=1 to !sb_k
      for !sb_l=1 to !sb_k
        !sb_hat=!sb_hat+w(!sb_i,!sb_j)*gi(!sb_j,!sb_l)*w(!sb_i,!sb_l)
      next
    next
    !sb_d=1-!sb_hat
    if !sb_d<0.05 then
      !sb_d=0.05
    endif
    innov(!sb_i)=e(!sb_i)/@sqrt(!sb_d)
  next
  vector ys
  matrix ws
  vector vs
  scalar nds
  scalar f1
  scalar c1
  scalar s1
  scalar ccut
  matrix path1
  vector(reps) fdraw
  vector(reps) cdraw
  vector(reps) sdraw
  !sb_tf=0
  !sb_tc=0
  !sb_ts=0
  for !sb_b=1 to reps
    if !sb_b=1 or (!sb_b-1)=25*@floor((!sb_b-1)/25) then
      call progress("Recursive stability bootstrap",!sb_b-1,reps,model_id)
    endif
    for !sb_i=1 to !sb_n
      eps(!sb_i)=innov(!sb_i)*(2*(@runif(0,1)>=0.5)-1)
    next
    call recurse_y(y,x,ord,det,hold,w,nd,b,eps,ys)
    call refresh_ecm(ys,ord,hold,w,nd,ws,vs)
      nds=nd
    call supf(ws,vs,f1,ccut)
    call recursive_paths(ws,vs,first,c1,s1,path1)
    fdraw(!sb_b)=f1
    cdraw(!sb_b)=c1
    sdraw(!sb_b)=s1
    if cut>=0 then
      !sb_tf=!sb_tf+(f1>=fs)
    endif
    !sb_tc=!sb_tc+(c1>=cm)
    !sb_ts=!sb_ts+(s1>=qm)
  next
  answer(3)=(1+!sb_tf)/(reps+1)
  if cut<0 then
    answer(3)=NA
  endif
  answer(5)=(1+!sb_tc)/(reps+1)
  answer(7)=(1+!sb_ts)/(reps+1)
  scalar quant
  call quantile7(cdraw,0.95,quant)
  answer(9)=quant
  call quantile7(sdraw,0.95,quant)
  answer(10)=quant
endsub

subroutine local adf_test(vector z, scalar trendflag, scalar maxlag, scalar bestlag, scalar adfstat, scalar adfp)
  !ad_n=@rows(z)
  !ad_best=1e100
  matrix w
  vector v
  vector b
  vector e
  matrix gi
  scalar ss
  scalar rk
  for !ad_l=0 to maxlag
    !ad_nr=!ad_n-maxlag-1
    w=@ones(!ad_nr,2+trendflag+!ad_l)
    v=@zeros(!ad_nr)
    for !ad_i=1 to !ad_nr
      !ad_t=maxlag+1+!ad_i
      v(!ad_i)=z(!ad_t)-z(!ad_t-1)
      w(!ad_i,2)=z(!ad_t-1)
      if trendflag=1 then
        w(!ad_i,3)=!ad_i
      endif
      if !ad_l>0 then
        for !ad_j=1 to !ad_l
          w(!ad_i,2+trendflag+!ad_j)=z(!ad_t-!ad_j)-z(!ad_t-!ad_j-1)
        next
      endif
    next
    call olsfit(w,v,b,e,gi,ss,rk)
    !ad_bic=!ad_nr*(log(2*3.141592653589793)+1+log(ss/!ad_nr))+@columns(w)*log(!ad_nr)
    if !ad_bic<!ad_best then
      !ad_best=!ad_bic
      bestlag=!ad_l
    endif
  next
  !ad_nr=!ad_n-bestlag-1
  w=@ones(!ad_nr,2+trendflag+bestlag)
  v=@zeros(!ad_nr)
  for !ad_i=1 to !ad_nr
    !ad_t=bestlag+1+!ad_i
    v(!ad_i)=z(!ad_t)-z(!ad_t-1)
    w(!ad_i,2)=z(!ad_t-1)
    if trendflag=1 then
      w(!ad_i,3)=!ad_i
    endif
    if bestlag>0 then
      for !ad_j=1 to bestlag
        w(!ad_i,2+trendflag+!ad_j)=z(!ad_t-!ad_j)-z(!ad_t-!ad_j-1)
      next
    endif
  next
  call olsfit(w,v,b,e,gi,ss,rk)
  adfstat=b(2)/@sqrt(gi(2,2)*ss/(!ad_nr-@columns(w)))
  ' MacKinnon (1994) N=1 response surfaces, c and ct.
  if trendflag=0 then
    !ad_min=-18.83
    !ad_max=2.74
    if adfstat<=-1.61 then
      !ad_poly=2.1659+1.4412*adfstat+0.038269*adfstat^2
    else
      !ad_poly=1.7339+0.93202*adfstat-0.12745*adfstat^2-0.010368*adfstat^3
    endif
  else
    !ad_min=-16.18
    !ad_max=0.7
    if adfstat<=-2.89 then
      !ad_poly=3.2512+1.6047*adfstat+0.049588*adfstat^2
    else
      !ad_poly=2.5261+0.61654*adfstat-0.37956*adfstat^2-0.060285*adfstat^3
    endif
  endif
  adfp=@cnorm(!ad_poly)
  if adfstat<!ad_min then
    adfp=0
  endif
  if adfstat>!ad_max then
    adfp=1
  endif
endsub

subroutine local kpss_test(vector z, scalar trendflag, scalar stat, scalar prob, scalar bandwidth)
  !kp_n=@rows(z)
  matrix w=@ones(!kp_n,1+trendflag)
  if trendflag=1 then
    for !kp_i=1 to !kp_n
      w(!kp_i,2)=!kp_i
    next
  endif
  vector b
  vector e
  matrix gi
  scalar ss
  scalar rk
  call olsfit(w,z,b,e,gi,ss,rk)
  !kp_covlags=@floor(!kp_n^(2/9))
  !kp_s0=ss/!kp_n
  !kp_s1=0
  for !kp_l=1 to !kp_covlags
    !kp_prod=0
    for !kp_i=!kp_l+1 to !kp_n
      !kp_prod=!kp_prod+e(!kp_i)*e(!kp_i-!kp_l)
    next
    !kp_prod=2*!kp_prod/!kp_n
    !kp_s0=!kp_s0+!kp_prod
    !kp_s1=!kp_s1+!kp_l*!kp_prod
  next
  bandwidth=@floor(1.1447*((!kp_s1/!kp_s0)^2)^(1/3)*!kp_n^(1/3))
  if bandwidth>!kp_n-1 then
    bandwidth=!kp_n-1
  endif
  !kp_var=ss
  if bandwidth>0 then
    for !kp_l=1 to bandwidth
      !kp_prod=0
      for !kp_i=!kp_l+1 to !kp_n
        !kp_prod=!kp_prod+e(!kp_i)*e(!kp_i-!kp_l)
      next
      !kp_var=!kp_var+2*(1-!kp_l/(bandwidth+1))*!kp_prod
    next
  endif
  !kp_var=!kp_var/!kp_n
  !kp_cum=0
  !kp_eta=0
  for !kp_i=1 to !kp_n
    !kp_cum=!kp_cum+e(!kp_i)
    !kp_eta=!kp_eta+!kp_cum^2
  next
  stat=!kp_eta/!kp_n^2/!kp_var
  vector(4) cv
  vector(4) pv
  pv(1)=0.10
  pv(2)=0.05
  pv(3)=0.025
  pv(4)=0.01
  if trendflag=0 then
    cv(1)=0.347
    cv(2)=0.463
    cv(3)=0.574
    cv(4)=0.739
  else
    cv(1)=0.119
    cv(2)=0.146
    cv(3)=0.176
    cv(4)=0.216
  endif
  prob=0.10
  if stat>=cv(4) then
    prob=0.01
  else
    for !kp_j=1 to 3
      if stat>=cv(!kp_j) and stat<cv(!kp_j+1) then
        prob=pv(!kp_j)+(stat-cv(!kp_j))/(cv(!kp_j+1)-cv(!kp_j))*(pv(!kp_j+1)-pv(!kp_j))
      endif
    next
  endif
endsub

subroutine local za_test(vector z, matrix lookup, scalar stat, scalar prob, scalar lag, scalar breakyear)
  scalar ast
  scalar ap
  call adf_test(z,1,4,lag,ast,ap)
  !za_n=@rows(z)
  !za_nr=!za_n-lag-1
  matrix w=@zeros(!za_nr,5+lag)
  vector v=@zeros(!za_nr)
  vector b
  vector e
  matrix gi
  scalar ss
  scalar rk
  stat=1e100
  for !za_bp=@floor(0.15*!za_n)+1 to !za_n-@floor(0.15*!za_n)
    for !za_i=1 to !za_nr
      !za_t=lag+1+!za_i
      v(!za_i)=z(!za_t)-z(!za_t-1)
      w(!za_i,1)=1
      w(!za_i,2)=(!za_t>!za_bp)
      w(!za_i,3)=!za_t+1
      w(!za_i,4)=(!za_t>!za_bp)*(!za_t-!za_bp+1)
      w(!za_i,5)=z(!za_t-1)
      if lag>0 then
        for !za_j=1 to lag
          w(!za_i,5+!za_j)=z(!za_t-!za_j)-z(!za_t-!za_j-1)
        next
      endif
    next
    call olsfit(w,v,b,e,gi,ss,rk)
    !za_test=b(5)/@sqrt(gi(5,5)*ss/(!za_nr-@columns(w)))
    if !za_test<stat then
      stat=!za_test
      breakyear=1964+!za_bp
    endif
  next
  prob=lookup(1,1)/100
  if stat>=lookup(@rows(lookup),2) then
    prob=lookup(@rows(lookup),1)/100
  else
    for !za_j=1 to @rows(lookup)-1
      if stat>=lookup(!za_j,2) and stat<lookup(!za_j+1,2) then
        prob=(lookup(!za_j,1)+(stat-lookup(!za_j,2))/(lookup(!za_j+1,2)-lookup(!za_j,2))*(lookup(!za_j+1,1)-lookup(!za_j,1)))/100
      endif
    next
  endif
endsub

subroutine local pss_reference(scalar n, scalar reps, matrix distributions)
  ' I(1) upper-reference distributions for Cases 3,4,5; exact n, five x's.
  distributions=@zeros(reps,3)
  matrix innovations
  matrix walk
  matrix w
  vector v
  matrix wr
  vector b
  vector e
  matrix gi
  scalar ss
  scalar rk
  scalar sr
  matrix accumulator=@zeros(n+1,n+1)
  for !ps_i=1 to n+1
    for !ps_j=1 to !ps_i
      accumulator(!ps_i,!ps_j)=1
    next
  next
  vector timeindex=@zeros(n)
  for !ps_i=1 to n
    timeindex(!ps_i)=!ps_i-1
  next
  matrix levelblock
  for !ps_case=3 to 5
    rndseed 46353
    for !ps_b=1 to reps
      if !ps_b=1 or (!ps_b-1)=500*@floor((!ps_b-1)/500) then
        %ps_msg="PSS simulation case="+@str(!ps_case)
        call progress(%ps_msg,!ps_b-1,reps,0)
      endif
      innovations=@mnrnd(n+1,6)
      walk=accumulator*innovations
      !ps_k=7+(!ps_case>3)
      w=@ones(n,!ps_k)
      v=@subextract(walk,2,1,n+1,1)-@subextract(walk,1,1,n,1)
      levelblock=@subextract(walk,1,1,n,6)
      matplace(w,levelblock,1,1)
      if !ps_case>3 then
        colplace(w,timeindex,7)
      endif
      call olsfit(w,v,b,e,gi,ss,rk)
      !ps_sse=ss
      wr=@ones(n,1)
      if !ps_case=5 then
        wr=@ones(n,2)
        colplace(wr,timeindex,1)
      endif
      call olsfit(wr,v,b,e,gi,sr,rk)
      !ps_r=!ps_k-@columns(wr)
      distributions(!ps_b,!ps_case-2)=((sr-!ps_sse)/!ps_r)/(!ps_sse/(n-!ps_k))
    next
  next
endsub

subroutine local influence_delete(matrix w, vector v, scalar nd, scalar hold, matrix results)
  !id_n=@rows(w)
  !id_k=@columns(w)
  results=@zeros(2*!id_n-4,6)
  matrix wd
  vector vd
  vector b
  vector e
  matrix gi
  scalar ss
  scalar rk
  !id_row=0
  for !id_mode=1 to 2
    !id_width=1+4*(!id_mode=2)
    for !id_start=1 to !id_n-!id_width+1
      wd=@zeros(!id_n-!id_width,!id_k)
      vd=@zeros(!id_n-!id_width)
      !id_r=0
      for !id_i=1 to !id_n
        if !id_i<!id_start or !id_i>=!id_start+!id_width then
          !id_r=!id_r+1
          vd(!id_r)=v(!id_i)
          for !id_j=1 to !id_k
            wd(!id_r,!id_j)=w(!id_i,!id_j)
          next
        endif
      next
      call olsfit(wd,vd,b,e,gi,ss,rk)
      !id_row=!id_row+1
      results(!id_row,1)=!id_width
      results(!id_row,2)=1964+hold+!id_start
      results(!id_row,3)=1963+hold+!id_start+!id_width
      results(!id_row,4)=b(nd+6)
      !id_se=@sqrt(gi(nd+6,nd+6)*ss/(@rows(wd)-!id_k))
      results(!id_row,5)=2*(1-@ctdist(@abs(b(nd+6)/!id_se),@rows(wd)-!id_k))
      results(!id_row,6)=b(nd+1)
    next
  next
endsub


subroutine local refresh_ecm(vector ys, vector ord, scalar hold, matrix template, scalar nd, matrix ws, vector vs)
  !rf_T=@rows(ys)
  ws=template
  vector piece=@subextract(ys,hold,1,!rf_T-1,1)
  vs=@subextract(ys,hold+1,1,!rf_T,1)-piece
  colplace(ws,piece,nd+1)
  if ord(1)>1 then
    for !rf_l=1 to ord(1)-1
      piece=@subextract(ys,hold+1-!rf_l,1,!rf_T-!rf_l,1)-@subextract(ys,hold-!rf_l,1,!rf_T-!rf_l-1,1)
      colplace(ws,piece,nd+6+!rf_l)
    next
  endif
endsub

subroutine local progress(string %pg_stage, scalar current, scalar total, scalar model_id)
  table(8,2) heartbeat
  heartbeat(1,1)="Build"
  heartbeat(1,2)="2026-09-26-r6"
  heartbeat(2,1)="Stage"
  heartbeat(2,2)=%pg_stage
  heartbeat(3,1)="Model ID (0 during search)"
  heartbeat(3,2)=model_id
  heartbeat(4,1)="Completed within stage"
  heartbeat(4,2)=current
  heartbeat(5,1)="Total within stage"
  heartbeat(5,2)=total
  heartbeat(6,1)="Completed"
  heartbeat(6,2)=current
  heartbeat(7,1)="Status"
  heartbeat(7,2)="RUNNING; only LCF_STATUS.csv SUCCESS marks completion"
  if %pg_stage="All computations and exports completed" then
    heartbeat(7,2)="SUCCESS"
  endif
  heartbeat.save(t=csv) LCF_PROGRESS.csv
  %pg_message=%pg_stage+" | model "+@str(model_id)+" | "+@str(current)+" / "+@str(total)
  statusline %pg_message
  logmsg %pg_message
endsub

' ======================== STUDY EXECUTION ========================
' This is a full production run. No external process or saved estimates are used.
' Bootstrap streams are EViews-native; Monte Carlo p-values need not be bitwise
' identical to runs in other software. Seeds and replication counts are saved.
mode quiet
logmode l e -p -s
!run_model=0
!B=1999
!BS=999
!NPSS=30000
%outputdir=@runpath
if @lower(@right(%outputdir,4))=".prg" then
  for !pathi=@len(%outputdir) to 1 step -1
    if @mid(%outputdir,!pathi,1)="\" then
      %outputdir=@left(%outputdir,!pathi)
      exitloop
    endif
  next
endif
if @right(%outputdir,1)<>"\" then
  %outputdir=%outputdir+"\"
endif
cd %outputdir
wfcreate(wf=LCF_NATIVE_RESULTS,page=annual) a 1965 2024
spool results_spool
table run_status
run_status(1,1)="RUNNING: do not interpret partially written outputs as a completed run."
run_status.save(t=csv) LCF_STATUS.csv
call progress("OLS compatibility preflight",0,1,!run_model)
' Small deterministic checks before expensive work; no study estimates embedded.
matrix pf_x=@ones(60,2)
vector pf_y=@zeros(60)
for !pf_i=1 to 60
  pf_x(!pf_i,2)=!pf_i/60
  pf_y(!pf_i)=2+3*pf_x(!pf_i,2)
next
vector pf_b
vector pf_e
matrix pf_g
scalar pf_ss
scalar pf_rank
call olsfit(pf_x,pf_y,pf_b,pf_e,pf_g,pf_ss,pf_rank)
if @abs(pf_b(1)-2)>1e-8 or @abs(pf_b(2)-3)>1e-8 or pf_ss>1e-12 or pf_rank<>2 then
  seterr "OLS preflight failed for two-column design."
endif
pf_x=@ones(60,1)
call olsfit(pf_x,pf_y,pf_b,pf_e,pf_g,pf_ss,pf_rank)
if @abs(pf_b(1)-@mean(pf_y))>1e-8 or @abs(pf_g(1,1)-1/60)>1e-10 then
  seterr "OLS preflight failed for single-column design."
endif
call progress("OLS preflight passed; loading embedded data",1,1,!run_model)
' Embedded, unchanged DATA values; source workbook SHA256:
' d810b2a5d5f4161abb88a39cf982e3383ccac19354ea43d5329a2ef35af95edb
smpl @all
series lcf=NA
smpl 1965 1965
lcf=1.3618194779999999
smpl 1966 1966
lcf=1.3318786739999999
smpl 1967 1967
lcf=1.363930219
smpl 1968 1968
lcf=1.2705773330000001
smpl 1969 1969
lcf=1.306804152
smpl 1970 1970
lcf=1.2642488759999999
smpl 1971 1971
lcf=1.18830535
smpl 1972 1972
lcf=1.263731648
smpl 1973 1973
lcf=1.286353753
smpl 1974 1974
lcf=1.299791658
smpl 1975 1975
lcf=0.96290713699999997
smpl 1976 1976
lcf=0.85539219200000005
smpl 1977 1977
lcf=1.1157378309999999
smpl 1978 1978
lcf=0.89937595999999997
smpl 1979 1979
lcf=0.7511563
smpl 1980 1980
lcf=0.67076232700000005
smpl 1981 1981
lcf=0.60332790199999997
smpl 1982 1982
lcf=0.61468066099999996
smpl 1983 1983
lcf=0.55448213700000004
smpl 1984 1984
lcf=0.52923203699999999
smpl 1985 1985
lcf=0.51256110700000002
smpl 1986 1986
lcf=0.54223068200000002
smpl 1987 1987
lcf=0.52000349599999995
smpl 1988 1988
lcf=0.52818181500000005
smpl 1989 1989
lcf=0.45379879499999998
smpl 1990 1990
lcf=0.478049893
smpl 1991 1991
lcf=0.47620567800000002
smpl 1992 1992
lcf=0.48431403699999998
smpl 1993 1993
lcf=0.47432438199999999
smpl 1994 1994
lcf=0.43917411299999998
smpl 1995 1995
lcf=0.42016826099999999
smpl 1996 1996
lcf=0.40031144000000002
smpl 1997 1997
lcf=0.40170127700000002
smpl 1998 1998
lcf=0.43458499299999998
smpl 1999 1999
lcf=0.36519322999999998
smpl 2000 2000
lcf=0.33355178000000002
smpl 2001 2001
lcf=0.32603163600000001
smpl 2002 2002
lcf=0.351172225
smpl 2003 2003
lcf=0.34642305400000001
smpl 2004 2004
lcf=0.32337274199999999
smpl 2005 2005
lcf=0.300867675
smpl 2006 2006
lcf=0.30333681099999998
smpl 2007 2007
lcf=0.29105588599999999
smpl 2008 2008
lcf=0.237044799
smpl 2009 2009
lcf=0.24522074799999999
smpl 2010 2010
lcf=0.24975813699999999
smpl 2011 2011
lcf=0.22996367700000001
smpl 2012 2012
lcf=0.22240511099999999
smpl 2013 2013
lcf=0.22448184199999999
smpl 2014 2014
lcf=0.227981656
smpl 2015 2015
lcf=0.25106678999999998
smpl 2016 2016
lcf=0.27303724800000001
smpl 2017 2017
lcf=0.244210028
smpl 2018 2018
lcf=0.25678617999999998
smpl 2019 2019
lcf=0.234267634
smpl 2020 2020
lcf=0.23050975800000001
smpl 2021 2021
lcf=0.21609605900000001
smpl 2022 2022
lcf=0.19848539800000001
smpl 2023 2023
lcf=0.190972635
smpl 2024 2024
lcf=0.188202967
smpl @all
series gdp_pc_real=NA
smpl 1965 1965
gdp_pc_real=3255.501358
smpl 1966 1966
gdp_pc_real=3518.6162180000001
smpl 1967 1967
gdp_pc_real=3797.076967
smpl 1968 1968
gdp_pc_real=4215.9367320000001
smpl 1969 1969
gdp_pc_real=4726.0925219999999
smpl 1970 1970
gdp_pc_real=5085.7568170000004
smpl 1971 1971
gdp_pc_real=5612.3630149999999
smpl 1972 1972
gdp_pc_real=6240.7007739999999
smpl 1973 1973
gdp_pc_real=6511.6545429999996
smpl 1974 1974
gdp_pc_real=6688.445213
smpl 1975 1975
gdp_pc_real=6468.1223950000003
smpl 1976 1976
gdp_pc_real=7422.1286090000003
smpl 1977 1977
gdp_pc_real=6996.1843140000001
smpl 1978 1978
gdp_pc_real=5901.8129600000002
smpl 1979 1979
gdp_pc_real=5016.9148679999998
smpl 1980 1980
gdp_pc_real=3793.3542969999999
smpl 1981 1981
gdp_pc_real=3394.343582
smpl 1982 1982
gdp_pc_real=3965.2755050000001
smpl 1983 1983
gdp_pc_real=4236.7168110000002
smpl 1984 1984
gdp_pc_real=3788.977629
smpl 1985 1985
gdp_pc_real=3720.5326230000001
smpl 1986 1986
gdp_pc_real=3234.8854310000002
smpl 1987 1987
gdp_pc_real=3115.6180319999999
smpl 1988 1988
gdp_pc_real=2833.187578
smpl 1989 1989
gdp_pc_real=2920.3693109999999
smpl 1990 1990
gdp_pc_real=3222.2010270000001
smpl 1991 1991
gdp_pc_real=3534.6539640000001
smpl 1992 1992
gdp_pc_real=3579.7322180000001
smpl 1993 1993
gdp_pc_real=3514.5196660000001
smpl 1994 1994
gdp_pc_real=3448.2595689999998
smpl 1995 1995
gdp_pc_real=3492.1720500000001
smpl 1996 1996
gdp_pc_real=3674.032686
smpl 1997 1997
gdp_pc_real=3679.275349
smpl 1998 1998
gdp_pc_real=3704.5656300000001
smpl 1999 1999
gdp_pc_real=3725.317458
smpl 2000 2000
gdp_pc_real=3885.3018339999999
smpl 2001 2001
gdp_pc_real=3917.2993550000001
smpl 2002 2002
gdp_pc_real=4198.9535610000003
smpl 2003 2003
gdp_pc_real=4511.4313400000001
smpl 2004 2004
gdp_pc_real=4607.4516910000002
smpl 2005 2005
gdp_pc_real=4650.3584149999997
smpl 2006 2006
gdp_pc_real=4778.7934580000001
smpl 2007 2007
gdp_pc_real=5084.7179599999999
smpl 2008 2008
gdp_pc_real=5035.9400850000002
smpl 2009 2009
gdp_pc_real=5023.8931439999997
smpl 2010 2010
gdp_pc_real=5249.061334
smpl 2011 2011
gdp_pc_real=5321.758707
smpl 2012 2012
gdp_pc_real=5058.6463780000004
smpl 2013 2013
gdp_pc_real=4916.9729189999998
smpl 2014 2014
gdp_pc_real=5093.2031969999998
smpl 2015 2015
gdp_pc_real=4952.7335549999998
smpl 2016 2016
gdp_pc_real=5312.6172500000002
smpl 2017 2017
gdp_pc_real=5395.2398400000002
smpl 2018 2018
gdp_pc_real=5127.2249540000003
smpl 2019 2019
gdp_pc_real=4952.4566910000003
smpl 2020 2020
gdp_pc_real=5132.8243080000002
smpl 2021 2021
gdp_pc_real=5300.6227200000003
smpl 2022 2022
gdp_pc_real=5465.3144750000001
smpl 2023 2023
gdp_pc_real=5687.8439630000003
smpl 2024 2024
gdp_pc_real=5834.4429689999997
smpl @all
series gfcf_gdp=NA
smpl 1965 1965
gfcf_gdp=31.283483050000001
smpl 1966 1966
gfcf_gdp=29.288941439999999
smpl 1967 1967
gfcf_gdp=34.435068350000002
smpl 1968 1968
gfcf_gdp=34.20136557
smpl 1969 1969
gfcf_gdp=36.28996111
smpl 1970 1970
gfcf_gdp=36.81983941
smpl 1971 1971
gfcf_gdp=34.757050509999999
smpl 1972 1972
gfcf_gdp=36.913162489999998
smpl 1973 1973
gfcf_gdp=34.21119522
smpl 1974 1974
gfcf_gdp=27.965169970000002
smpl 1975 1975
gfcf_gdp=42.205130390000001
smpl 1976 1976
gfcf_gdp=49.316200969999997
smpl 1977 1977
gfcf_gdp=49.133524610000002
smpl 1978 1978
gfcf_gdp=47.257904879999998
smpl 1979 1979
gfcf_gdp=31.103984700000002
smpl 1980 1980
gfcf_gdp=36.922323800000001
smpl 1981 1981
gfcf_gdp=33.604468609999998
smpl 1982 1982
gfcf_gdp=30.948544770000002
smpl 1983 1983
gfcf_gdp=39.184908759999999
smpl 1984 1984
gfcf_gdp=39.435324569999999
smpl 1985 1985
gfcf_gdp=32.878526479999998
smpl 1986 1986
gfcf_gdp=29.969787140000001
smpl 1987 1987
gfcf_gdp=26.26023679
smpl 1988 1988
gfcf_gdp=24.59869986
smpl 1989 1989
gfcf_gdp=24.661248860000001
smpl 1990 1990
gfcf_gdp=26.98318321
smpl 1991 1991
gfcf_gdp=34.674870910000003
smpl 1992 1992
gfcf_gdp=33.016896410000001
smpl 1993 1993
gfcf_gdp=27.827960390000001
smpl 1994 1994
gfcf_gdp=25.964357530000001
smpl 1995 1995
gfcf_gdp=24.688894269999999
smpl 1996 1996
gfcf_gdp=31.464518909999999
smpl 1997 1997
gfcf_gdp=32.650233319999998
smpl 1998 1998
gfcf_gdp=32.638353340000002
smpl 1999 1999
gfcf_gdp=31.844273489999999
smpl 2000 2000
gfcf_gdp=31.305453920000001
smpl 2001 2001
gfcf_gdp=35.799397540000001
smpl 2002 2002
gfcf_gdp=33.245783260000003
smpl 2003 2003
gfcf_gdp=32.658082700000001
smpl 2004 2004
gfcf_gdp=31.23817549
smpl 2005 2005
gfcf_gdp=28.557168770000001
smpl 2006 2006
gfcf_gdp=27.247096150000001
smpl 2007 2007
gfcf_gdp=28.342123999999998
smpl 2008 2008
gfcf_gdp=31.922693049999999
smpl 2009 2009
gfcf_gdp=30.85489639
smpl 2010 2010
gfcf_gdp=27.352833319999998
smpl 2011 2011
gfcf_gdp=29.05372873
smpl 2012 2012
gfcf_gdp=30.053965430000002
smpl 2013 2013
gfcf_gdp=27.44560637
smpl 2014 2014
gfcf_gdp=28.570813399999999
smpl 2015 2015
gfcf_gdp=25.365013900000001
smpl 2016 2016
gfcf_gdp=23.816045209999999
smpl 2017 2017
gfcf_gdp=23.97438189
smpl 2018 2018
gfcf_gdp=24.532164269999999
smpl 2019 2019
gfcf_gdp=25.24461307
smpl 2020 2020
gfcf_gdp=30.25955381
smpl 2021 2021
gfcf_gdp=28.62302764
smpl 2022 2022
gfcf_gdp=26.64418701
smpl 2023 2023
gfcf_gdp=27.416322189999999
smpl 2024 2024
gfcf_gdp=28.107160230000002
smpl @all
series energy_pc=NA
smpl 1965 1965
energy_pc=3885.8681999999999
smpl 1966 1966
energy_pc=4105.2969999999996
smpl 1967 1967
energy_pc=4380.7060000000001
smpl 1968 1968
energy_pc=4674.2924999999996
smpl 1969 1969
energy_pc=4978.3495999999996
smpl 1970 1970
energy_pc=5697.0537000000004
smpl 1971 1971
energy_pc=6056.8180000000002
smpl 1972 1972
energy_pc=6669.1356999999998
smpl 1973 1973
energy_pc=7740.9614000000001
smpl 1974 1974
energy_pc=8521.2739999999994
smpl 1975 1975
energy_pc=9633.732
smpl 1976 1976
energy_pc=10519.447
smpl 1977 1977
energy_pc=11547.531000000001
smpl 1978 1978
energy_pc=10926.727999999999
smpl 1979 1979
energy_pc=11398.028
smpl 1980 1980
energy_pc=10191.557000000001
smpl 1981 1981
energy_pc=9705.7240000000002
smpl 1982 1982
energy_pc=10389.114
smpl 1983 1983
energy_pc=11864.153
smpl 1984 1984
energy_pc=12348.089
smpl 1985 1985
energy_pc=12987.325999999999
smpl 1986 1986
energy_pc=11726.332
smpl 1987 1987
energy_pc=12128.986000000001
smpl 1988 1988
energy_pc=11976.522000000001
smpl 1989 1989
energy_pc=12714.513000000001
smpl 1990 1990
energy_pc=13896.569
smpl 1991 1991
energy_pc=14673.674999999999
smpl 1992 1992
energy_pc=15909.522000000001
smpl 1993 1993
energy_pc=14669.319
smpl 1994 1994
energy_pc=16396.866999999998
smpl 1995 1995
energy_pc=17007.226999999999
smpl 1996 1996
energy_pc=18302.474999999999
smpl 1997 1997
energy_pc=18638.291000000001
smpl 1998 1998
energy_pc=19149.243999999999
smpl 1999 1999
energy_pc=20346.25
smpl 2000 2000
energy_pc=21026.740000000002
smpl 2001 2001
energy_pc=22121.530999999999
smpl 2002 2002
energy_pc=23913.190999999999
smpl 2003 2003
energy_pc=24083.326000000001
smpl 2004 2004
energy_pc=26075.348000000002
smpl 2005 2005
energy_pc=27012.970000000001
smpl 2006 2006
energy_pc=28889.373
smpl 2007 2007
energy_pc=30213.291000000001
smpl 2008 2008
energy_pc=30767.291000000001
smpl 2009 2009
energy_pc=31118.98
smpl 2010 2010
energy_pc=31753.627
smpl 2011 2011
energy_pc=32584.848000000002
smpl 2012 2012
energy_pc=32428.09
smpl 2013 2013
energy_pc=33255.663999999997
smpl 2014 2014
energy_pc=34196.25
smpl 2015 2015
energy_pc=33441.343999999997
smpl 2016 2016
energy_pc=34845.239999999998
smpl 2017 2017
energy_pc=35561.832000000002
smpl 2018 2018
energy_pc=36251.887000000002
smpl 2019 2019
energy_pc=37875.233999999997
smpl 2020 2020
energy_pc=39001.836000000003
smpl 2021 2021
energy_pc=38159.269999999997
smpl 2022 2022
energy_pc=38985.188000000002
smpl 2023 2023
energy_pc=39688.279999999999
smpl 2024 2024
energy_pc=39203.610000000001
smpl @all
series sanction=NA
smpl 1965 1965
sanction=0
smpl 1966 1966
sanction=0
smpl 1967 1967
sanction=0
smpl 1968 1968
sanction=0
smpl 1969 1969
sanction=0
smpl 1970 1970
sanction=0
smpl 1971 1971
sanction=0
smpl 1972 1972
sanction=0
smpl 1973 1973
sanction=0
smpl 1974 1974
sanction=0
smpl 1975 1975
sanction=0
smpl 1976 1976
sanction=0
smpl 1977 1977
sanction=0
smpl 1978 1978
sanction=0
smpl 1979 1979
sanction=1
smpl 1980 1980
sanction=1
smpl 1981 1981
sanction=1
smpl 1982 1982
sanction=1
smpl 1983 1983
sanction=1
smpl 1984 1984
sanction=1
smpl 1985 1985
sanction=1
smpl 1986 1986
sanction=1
smpl 1987 1987
sanction=2
smpl 1988 1988
sanction=2
smpl 1989 1989
sanction=2
smpl 1990 1990
sanction=2
smpl 1991 1991
sanction=2
smpl 1992 1992
sanction=2
smpl 1993 1993
sanction=2
smpl 1994 1994
sanction=2
smpl 1995 1995
sanction=3
smpl 1996 1996
sanction=3
smpl 1997 1997
sanction=3
smpl 1998 1998
sanction=3
smpl 1999 1999
sanction=3
smpl 2000 2000
sanction=3
smpl 2001 2001
sanction=3
smpl 2002 2002
sanction=3
smpl 2003 2003
sanction=3
smpl 2004 2004
sanction=3
smpl 2005 2005
sanction=3
smpl 2006 2006
sanction=4
smpl 2007 2007
sanction=4
smpl 2008 2008
sanction=4
smpl 2009 2009
sanction=4
smpl 2010 2010
sanction=4
smpl 2011 2011
sanction=4
smpl 2012 2012
sanction=5
smpl 2013 2013
sanction=5
smpl 2014 2014
sanction=5
smpl 2015 2015
sanction=5
smpl 2016 2016
sanction=2
smpl 2017 2017
sanction=2
smpl 2018 2018
sanction=5
smpl 2019 2019
sanction=5
smpl 2020 2020
sanction=5
smpl 2021 2021
sanction=5
smpl 2022 2022
sanction=5
smpl 2023 2023
sanction=5
smpl 2024 2024
sanction=5

smpl @all
series lnlcf=log(lcf)
series lngdppc=log(gdp_pc_real)
series gfcf=gfcf_gdp
series lnenergy=log(energy_pc)
scalar sanction_mean=@mean(sanction)
scalar sanction_sd=@stdevp(sanction)
series si_z=(sanction-sanction_mean)/sanction_sd
series interaction=si_z*lnenergy
series trend=@trend/60
series post1979=(@trend+1965>=1979)
group gx lngdppc gfcf lnenergy si_z interaction
group gyx lnlcf lngdppc gfcf lnenergy si_z interaction
matrix xbase=@convert(gx)
matrix allseries=@convert(gyx)
vector y=@columnextract(allseries,1)
vector rawsan=@convert(sanction)
vector ord
vector obic
vector oaicc
matrix candidates
matrix mgrid=@zeros(74,12)
svector(74) model_names
!row=0
!primary=0
!trendmodel=0
%varnames="lnlcf lngdppc gfcf lnenergy si_z interaction"
%detnames="c ct cb ctb"
%critnames="BIC AICc"

' Integration tests: BIC ADF, matching Hobijn automatic KPSS, ct ZA.
matrix(48,2) za_lookup
za_lookup(1,1)=0.001
za_lookup(1,2)=-38.177999999999997
za_lookup(2,1)=0.10000000000000001
za_lookup(2,2)=-6.4310700000000001
za_lookup(3,1)=0.20000000000000001
za_lookup(3,2)=-6.0727900000000004
za_lookup(4,1)=0.29999999999999999
za_lookup(4,2)=-5.9549599999999998
za_lookup(5,1)=0.40000000000000002
za_lookup(5,2)=-5.8625400000000001
za_lookup(6,1)=0.5
za_lookup(6,2)=-5.77081
za_lookup(7,1)=0.59999999999999998
za_lookup(7,2)=-5.7254100000000001
za_lookup(8,1)=0.69999999999999996
za_lookup(8,2)=-5.6840599999999997
za_lookup(9,1)=0.80000000000000004
za_lookup(9,2)=-5.6516299999999999
za_lookup(10,1)=0.90000000000000002
za_lookup(10,2)=-5.60419
za_lookup(11,1)=1
za_lookup(11,2)=-5.5755600000000003
za_lookup(12,1)=2.5
za_lookup(12,2)=-5.29704
za_lookup(13,1)=5
za_lookup(13,2)=-5.0733199999999998
za_lookup(14,1)=7.5
za_lookup(14,2)=-4.9300300000000004
za_lookup(15,1)=10
za_lookup(15,2)=-4.8266799999999996
za_lookup(16,1)=12.5
za_lookup(16,2)=-4.7371100000000004
za_lookup(17,1)=15
za_lookup(17,2)=-4.6601999999999997
za_lookup(18,1)=17.5
za_lookup(18,2)=-4.5896999999999997
za_lookup(19,1)=20
za_lookup(19,2)=-4.5285500000000001
za_lookup(20,1)=22.5
za_lookup(20,2)=-4.4710000000000001
za_lookup(21,1)=25
za_lookup(21,2)=-4.4201100000000002
za_lookup(22,1)=27.5
za_lookup(22,2)=-4.3738700000000001
za_lookup(23,1)=30
za_lookup(23,2)=-4.3270499999999998
za_lookup(24,1)=32.5
za_lookup(24,2)=-4.2812599999999996
za_lookup(25,1)=35
za_lookup(25,2)=-4.2379300000000004
za_lookup(26,1)=37.5
za_lookup(26,2)=-4.1982200000000001
za_lookup(27,1)=40
za_lookup(27,2)=-4.1580000000000004
za_lookup(28,1)=42.5
za_lookup(28,2)=-4.1194600000000001
za_lookup(29,1)=45
za_lookup(29,2)=-4.0806399999999998
za_lookup(30,1)=47.5
za_lookup(30,2)=-4.0428600000000001
za_lookup(31,1)=50
za_lookup(31,2)=-4.0048899999999996
za_lookup(32,1)=52.5
za_lookup(32,2)=-3.9683700000000002
za_lookup(33,1)=55
za_lookup(33,2)=-3.9319999999999999
za_lookup(34,1)=57.5
za_lookup(34,2)=-3.8949600000000002
za_lookup(35,1)=60
za_lookup(35,2)=-3.8557700000000001
za_lookup(36,1)=65
za_lookup(36,2)=-3.7779500000000001
za_lookup(37,1)=70
za_lookup(37,2)=-3.69794
za_lookup(38,1)=75
za_lookup(38,2)=-3.6185200000000002
za_lookup(39,1)=80
za_lookup(39,2)=-3.5248499999999998
za_lookup(40,1)=85
za_lookup(40,2)=-3.4166500000000002
za_lookup(41,1)=90
za_lookup(41,2)=-3.2852700000000001
za_lookup(42,1)=92.5
za_lookup(42,2)=-3.1972399999999999
za_lookup(43,1)=95
za_lookup(43,2)=-3.0876899999999998
za_lookup(44,1)=96
za_lookup(44,2)=-3.0308799999999998
za_lookup(45,1)=97
za_lookup(45,2)=-2.9609100000000002
za_lookup(46,1)=98
za_lookup(46,2)=-2.85581
za_lookup(47,1)=99
za_lookup(47,2)=-2.7101500000000001
za_lookup(48,1)=99.900000000000006
za_lookup(48,2)=-2.2876699999999999

matrix units=@zeros(72,10)
matrix za_results=@zeros(6,5)
vector zz
vector workz
scalar alag
scalar astat
scalar ap
scalar kstat
scalar kp
scalar klag
scalar zstat
scalar zp
scalar zlag
scalar zbreak
!ur=0
for !uv=1 to 6
  call progress("Integration tests",!uv-1,6,!run_model)
  zz=@columnextract(allseries,!uv)
  for !ud=0 to 1
    workz=zz
    if !ud=1 then
      workz=@zeros(59)
      for !ui=1 to 59
        workz(!ui)=zz(!ui+1)-zz(!ui)
      next
    endif
    for !ut=0 to 1
      call kpss_test(workz,!ut,kstat,kp,klag)
      for !ul=2 to 6 step 2
        call adf_test(workz,!ut,!ul,alag,astat,ap)
        !ur=!ur+1
        units(!ur,1)=!uv
        units(!ur,2)=!ud
        units(!ur,3)=!ut
        units(!ur,4)=!ul
        units(!ur,5)=astat
        units(!ur,6)=ap
        units(!ur,7)=alag
        units(!ur,8)=kstat
        units(!ur,9)=kp
        units(!ur,10)=klag
      next
    next
  next
  call za_test(zz,za_lookup,zstat,zp,zlag,zbreak)
  za_results(!uv,1)=!uv
  za_results(!uv,2)=zstat
  za_results(!uv,3)=zp
  za_results(!uv,4)=zlag
  za_results(!uv,5)=zbreak
next
table unit_table
table za_table
call export_table(units,"variable_id difference trend maxlag ADF_stat ADF_p ADF_lag KPSS_stat KPSS_p_bounded KPSS_lag","LCF_unit_roots.csv",unit_table)
call export_table(za_results,"variable_id ZA_stat ZA_p ZA_lag break_year","LCF_zivot_andrews.csv",za_table)
results_spool.append unit_table za_table
wfsave LCF_NATIVE_RESULTS.wf1

' Exhaustive hierarchical lag search; every candidate uses holdback=4.
for !cap=1 to 4
  for !det=1 to 4
    call lagsearch(y,xbase,!cap,!det,4,0,obic,oaicc,candidates)
    %searchfile="LCF_search_L"+@str(!cap)+"_"+@word(%detnames,!det)+".csv"
    table search_export
    call export_table(candidates,"p q_gdp q_gfcf q_energy q_si q_interaction n k loglik BIC AICc",%searchfile,search_export)
    for !crit=1 to 2
      !row=!row+1
      ord=obic
      if !crit=2 then
        ord=oaicc
      endif
      mgrid(!row,1)=1
      mgrid(!row,2)=!cap
      mgrid(!row,3)=!det
      mgrid(!row,4)=!crit
      mgrid(!row,5)=4
      for !j=1 to 6
        mgrid(!row,!j+5)=ord(!j)
      next
      model_names(!row)=@word(%critnames,!crit)+"_L"+@str(!cap)+"_"+@word(%detnames,!det)
      if !cap=4 and !det=1 and !crit=1 then
        !primary=!row
      endif
      if !cap=4 and !det=2 and !crit=1 then
        !trendmodel=!row
      endif
    next
  next
next
' Start-year sensitivity: independently selected L=3, intercept, holdback=3.
call lagsearch(y,xbase,3,1,3,0,obic,oaicc,candidates)
table start_search
call export_table(candidates,"p q_gdp q_gfcf q_energy q_si q_interaction n k loglik BIC AICc","LCF_search_start1968.csv",start_search)
for !crit=1 to 2
  !row=!row+1
  ord=obic
  if !crit=2 then
    ord=oaicc
  endif
  mgrid(!row,1)=2
  mgrid(!row,2)=3
  mgrid(!row,3)=1
  mgrid(!row,4)=!crit
  mgrid(!row,5)=3
  for !j=1 to 6
    mgrid(!row,!j+5)=ord(!j)
  next
  model_names(!row)=@word(%critnames,!crit)+"_start1968"
next
' All-q>=1 sensitivity retains the same grid and common sample.
for !cap=1 to 4
  for !det=1 to 4
    call lagsearch(y,xbase,!cap,!det,4,1,obic,oaicc,candidates)
    %searchfile="LCF_search_qmin1_L"+@str(!cap)+"_"+@word(%detnames,!det)+".csv"
    table qsearch_export
    call export_table(candidates,"p q_gdp q_gfcf q_energy q_si q_interaction n k loglik BIC AICc",%searchfile,qsearch_export)
    for !crit=1 to 2
      !row=!row+1
      ord=obic
      if !crit=2 then
        ord=oaicc
      endif
      mgrid(!row,1)=3
      mgrid(!row,2)=!cap
      mgrid(!row,3)=!det
      mgrid(!row,4)=!crit
      mgrid(!row,5)=4
      for !j=1 to 6
        mgrid(!row,!j+5)=ord(!j)
      next
      model_names(!row)="qmin1_"+@word(%critnames,!crit)+"_L"+@str(!cap)+"_"+@word(%detnames,!det)
    next
  next
next
' Ordinal recoding: restandardization is always over all 60 observations.
vector sc
vector zs
matrix xx
%codings="original sqrt square log1p"
for !code=0 to 3
  sc=rawsan
  if !code=1 then
    sc=@sqrt(rawsan)
  endif
  if !code=2 then
    sc=@epow(rawsan,2)
  endif
  if !code=3 then
    sc=@log(rawsan+1)
  endif
  zs=(sc-@mean(sc))/@sqrt(@mean(@epow(sc-@mean(sc),2)))
  xx=xbase
  for !i=1 to 60
    xx(!i,4)=zs(!i)
    xx(!i,5)=zs(!i)*xx(!i,3)
  next
  for !det=1 to 2
    call lagsearch(y,xx,3,!det,4,0,obic,oaicc,candidates)
    !row=!row+1
    mgrid(!row,1)=4
    mgrid(!row,2)=3
    mgrid(!row,3)=!det
    mgrid(!row,4)=1
    mgrid(!row,5)=4
    mgrid(!row,12)=!code
    for !j=1 to 6
      mgrid(!row,!j+5)=obic(!j)
    next
    model_names(!row)="coding_"+@word(%codings,!code+1)+"_"+@word(%detnames,!det)
  next
next
table model_index
call export_table(mgrid,"family maxlag deterministic criterion holdback p q_gdp q_gfcf q_energy q_si q_interaction coding","LCF_model_index.csv",model_index)
model_index(1,13)="model_id"
model_index(1,14)="model_name"
for !i=1 to 74
  model_index(!i+1,13)=!i
  model_index(!i+1,14)=model_names(!i)
next
model_index.save(t=csv) LCF_model_index.csv
results_spool.append model_index
wfsave LCF_NATIVE_RESULTS.wf1

' Selected models: native equation objects plus independently assembled exact ECM.
matrix summary=@zeros(74,28)
matrix diagall=@zeros(74,27)
matrix coefall=@zeros(2500,8)
matrix lrall=@zeros(1200,8)
!coefrow=0
!lrrow=0
matrix w
matrix a
vector v
vector ya
vector b
vector e
vector ba
vector ea
matrix gi
matrix ga
matrix cv
matrix hac
vector lev
vector cook
scalar ss
scalar sa
scalar rk
scalar rka
scalar nd
scalar nda
scalar armax
matrix boot
matrix lr
vector dg
vector phi
!seed_counter=0
for !m=1 to 74
  !run_model=!m
  call progress("Selected model started",!m-1,74,!run_model)
  !det=mgrid(!m,3)
  !hold=mgrid(!m,5)
  !family=mgrid(!m,1)
  !code=mgrid(!m,12)
  ord=@transpose(@subextract(mgrid,!m,6,!m,11))
  !reps=!B
  if !family=4 then
    !reps=!BS
  endif
  sc=rawsan
  if !code=1 then
    sc=@sqrt(rawsan)
  endif
  if !code=2 then
    sc=@epow(rawsan,2)
  endif
  if !code=3 then
    sc=@log(rawsan+1)
  endif
  zs=(sc-@mean(sc))/@sqrt(@mean(@epow(sc-@mean(sc),2)))
  xx=xbase
  for !i=1 to 60
    xx(!i,4)=zs(!i)
    xx(!i,5)=zs(!i)*xx(!i,3)
  next
  call design(y,xx,ord,!det,!hold,1,w,v,nd)
  call design(y,xx,ord,!det,!hold,0,a,ya,nda)
  call olsfit(w,v,b,e,gi,ss,rk)
  call olsfit(a,ya,ba,ea,ga,sa,rka)
  !equiv=@max(@abs(e-ea))
  if !equiv>1e-8 then
    seterr "Exact ARDL/ECM residual equality failed."
  endif
  !n=@rows(w)
  !k=@columns(w)
  !ll=-!n/2*(log(2*3.141592653589793)+1+log(sa/!n))
  summary(!m,1)=!m
  summary(!m,2)=!n
  summary(!m,3)=!k
  summary(!m,4)=-2*!ll+(!k+1)*log(!n)
  summary(!m,5)=-2*!ll+2*(!k+1)+2*(!k+1)*(!k+2)/(!n-!k-2)
  call covariances(w,e,gi,ss,cv,hac,lev,cook)
  summary(!m,6)=b(nd+1)
  summary(!m,7)=@sqrt(cv(nd+1,nd+1))
  summary(!m,8)=b(nd+6)
  summary(!m,9)=@sqrt(cv(nd+6,nd+6))
  summary(!m,10)=2*(1-@ctdist(@abs(b(nd+6)/summary(!m,9)),!n-!k))
  summary(!m,11)=2*(1-@ctdist(@abs(b(nd+6)/@sqrt(hac(nd+6,nd+6))),!n-!k))
  !cached=0
  if !m>1 then
    for !prev=1 to !m-1
      !same=(mgrid(!prev,3)=!det)*(mgrid(!prev,5)=!hold)*(mgrid(!prev,12)=!code)*(summary(!prev,28)=!reps)
      for !j=6 to 11
        !same=!same*(mgrid(!prev,!j)=mgrid(!m,!j))
      next
      if !same=1 then
        !cached=!prev
        exitloop
      endif
    next
  endif
  if !cached=0 then
    !seed=98237+!seed_counter
    !seed_counter=!seed_counter+1
    call bounds_boot(y,xx,ord,!det,!hold,!reps,!seed,0,!m,boot)
    for !j=1 to 3
      summary(!m,11+!j)=boot(!j,1)
      summary(!m,14+!j)=boot(!j,2)
      summary(!m,17+!j)=boot(!j,3)
      summary(!m,20+!j)=boot(!j,4)
    next
    summary(!m,27)=!seed
  else
    for !j=12 to 23
      summary(!m,!j)=summary(!cached,!j)
    next
    summary(!m,27)=summary(!cached,27)
  endif
  summary(!m,28)=!reps
  summary(!m,24)=(summary(!m,15)<0.05)*(summary(!m,16)<0.05)*(summary(!m,17)<0.05)*(b(nd+1)<0)
  phi=@subextract(ba,nd+1,1,nd+ord(1),1)
  call ar_modulus(phi,armax)
  summary(!m,25)=armax
  summary(!m,26)=!equiv
  call diagnostics(w,v,b,e,gi,ss,a,ya,ba,dg)
  diagall(!m,1)=!m
  for !j=1 to 24
    diagall(!m,!j+1)=dg(!j)
  next
  diagall(!m,26)=@max(lev)
  diagall(!m,27)=@max(cook)
  for !j=1 to !k
    !coefrow=!coefrow+1
    coefall(!coefrow,1)=!m
    coefall(!coefrow,2)=!j
    coefall(!coefrow,3)=b(!j)
    coefall(!coefrow,4)=@sqrt(cv(!j,!j))
    coefall(!coefrow,5)=2*(1-@ctdist(@abs(b(!j)/coefall(!coefrow,4)),!n-!k))
    coefall(!coefrow,6)=@sqrt(hac(!j,!j))
    coefall(!coefrow,7)=2*(1-@ctdist(@abs(b(!j)/coefall(!coefrow,6)),!n-!k))
    coefall(!coefrow,8)=!n-!k
  next
  if summary(!m,24)=1 then
    for !covtype=1 to 2
      if !covtype=1 then
        call long_run(b,cv,nd,!n-!k,lr)
      else
        call long_run(b,hac,nd,!n-!k,lr)
      endif
      for !term=1 to 8
        !lrrow=!lrrow+1
        lrall(!lrrow,1)=!m
        lrall(!lrrow,2)=!covtype
        lrall(!lrrow,3)=!term
        for !j=1 to 5
          lrall(!lrrow,!j+3)=lr(!term,!j)
        next
      next
    next
  endif
  ' Construct equation objects from the freshly selected order, never hard-coded.
  smpl @all
  %sivar="si_"+@str(!m)
  %intervar="iz_"+@str(!m)
  series {%sivar}
  mtos(zs,{%sivar})
  series {%intervar}={%sivar}*lnenergy
  %xnames="lngdppc gfcf lnenergy "+%sivar+" "+%intervar
  %dlist="c"
  if !det=2 or !det=4 then
    %dlist=%dlist+" trend"
  endif
  if !det=3 or !det=4 then
    %dlist=%dlist+" post1979"
  endif
  %alev=%dlist
  %ecml=%dlist+" lnlcf(-1)"
  for !l=1 to ord(1)
    %alev=%alev+" lnlcf(-"+@str(!l)+")"
  next
  for !j=1 to 5
    %xn=@word(%xnames,!j)
    for !l=0 to ord(!j+1)
      if !l=0 then
        %alev=%alev+" "+%xn
      else
        %alev=%alev+" "+%xn+"(-"+@str(!l)+")"
      endif
    next
    if ord(!j+1)=0 then
      %ecml=%ecml+" "+%xn
    else
      %ecml=%ecml+" "+%xn+"(-1)"
    endif
  next
  if ord(1)>1 then
    for !l=1 to ord(1)-1
      %ecml=%ecml+" d(lnlcf(-"+@str(!l)+"))"
    next
  endif
  for !j=1 to 5
    %xn=@word(%xnames,!j)
    if ord(!j+1)>0 then
      for !l=0 to ord(!j+1)-1
        if !l=0 then
          %ecml=%ecml+" d("+%xn+")"
        else
          %ecml=%ecml+" d("+%xn+"(-"+@str(!l)+"))"
        endif
      next
    endif
  next
  !firstyear=1965+!hold
  smpl {!firstyear} 2024
  %ae="ardl_"+@str(!m)
  %ee="ecm_"+@str(!m)
  equation {%ae}.ls lnlcf {%alev}
  equation {%ee}.ls d(lnlcf) {%ecml}
  for !j=1 to !k
    if @abs({%ee}.@coefs(!j)-b(!j))>1e-5 then
      seterr "Matrix/equation coefficient validation failed."
    endif
  next
  if !m=!primary or !m=!trendmodel then
    results_spool.append {%ae} {%ee}
  endif
  ' Only completed rows are exported in this explicitly partial checkpoint.
  table checkpoint_table
  matrix checkpoint_rows=@subextract(summary,1,1,!m,28)
  call export_table(checkpoint_rows,"model_id n k BIC AICc rho rho_se theta_inter theta_inter_se theta_inter_p HAC2_inter_p F t Fx pF pt pFx Fcrit5 tcrit5 Fxcrit5 mcseF mcset mcseFx level_support max_AR_modulus ECM_error bootstrap_seed B","LCF_models_PARTIAL.csv",checkpoint_table)
  call progress("Selected model completed",!m,74,!run_model)
  wfsave LCF_NATIVE_RESULTS.wf1
next
coefall=@subextract(coefall,1,1,!coefrow,8)
table main_table
table diagnostic_table
table coefficient_table
table longrun_table
call export_table(summary,"model_id n k BIC AICc rho rho_se theta_inter theta_inter_se theta_inter_p HAC2_inter_p F t Fx pF pt pFx Fcrit5 tcrit5 Fxcrit5 mcseF mcset mcseFx level_support max_AR_modulus ECM_error bootstrap_seed B","LCF_models.csv",main_table)
call export_table(diagall,"model_id BG1_p BG2_p BG4_p BG1_F_p BG2_F_p BG4_F_p LB1_p LB2_p LB4_p ARCH1_p ARCH2_p ARCH4_p BP_p JB_p DW RESET_delta_p RESET_level_p OLS_CUSUM OLS_CUSUM_p White_rank White_residual_df White_p condition_number residual_df max_leverage max_Cook","LCF_diagnostics.csv",diagnostic_table)
call export_table(coefall,"model_id coefficient_position estimate OLS_se OLS_p HAC2_se HAC2_p residual_df","LCF_ECM_coefficients.csv",coefficient_table)
if !lrrow>0 then
  lrall=@subextract(lrall,1,1,!lrrow,8)
  call export_table(lrall,"model_id covariance_id term_id estimate se lower95 upper95 delta_p","LCF_conditional_longrun.csv",longrun_table)
else
  longrun_table(1,1)="No model passed all three tests with negative adjustment. No long-run ratios reported."
  longrun_table.save(t=csv) LCF_conditional_longrun.csv
endif
results_spool.append main_table diagnostic_table coefficient_table longrun_table

' Conventional upper PSS reference: only all-q>=1 and no-shift models.
!run_model=0
call progress("Starting PSS reference simulations",0,3*!NPSS,!run_model)
matrix pss_draws
call pss_reference(56,!NPSS,pss_draws)
matrix pssout=@zeros(24,5)
!psrow=0
vector pdist
scalar pcrit
vector bl
matrix vl
matrix wal
for !m=35 to 66
  !det=mgrid(!m,3)
  if !det=1 or !det=2 then
    ord=@transpose(@subextract(mgrid,!m,6,!m,11))
    call design(y,xbase,ord,!det,4,1,w,v,nd)
    call olsfit(w,v,b,e,gi,ss,rk)
    !n=@rows(w)
    !k=@columns(w)
    !caselo=3
    !casehi=3
    if !det=2 then
      !caselo=4
      !casehi=5
    endif
    for !case=!caselo to !casehi
      !obsF=summary(!m,12)
      if !case=4 then
        bl=@subextract(b,2,1,nd+6,1)
        vl=@subextract(gi,2,2,nd+6,nd+6)*(ss/(!n-!k))
        wal=@transpose(bl)*@inverse(vl)*bl
        !obsF=wal(1,1)/7
      endif
      pdist=@columnextract(pss_draws,!case-2)
      !tail=0
      for !r=1 to !NPSS
        !tail=!tail+(pdist(!r)>!obsF)
      next
      call quantile7(pdist,0.95,pcrit)
      !psrow=!psrow+1
      pssout(!psrow,1)=!m
      pssout(!psrow,2)=!case
      pssout(!psrow,3)=!obsF
      pssout(!psrow,4)=!tail/!NPSS
      pssout(!psrow,5)=pcrit
    next
  endif
next
table pss_table
call export_table(pssout,"model_id case F upper_p upper_critical5","LCF_PSS_benchmarks.csv",pss_table)
results_spool.append pss_table
wfsave LCF_NATIVE_RESULTS.wf1

' IID sensitivity, stability, influence: primary and BIC trend counterpart.
matrix iidout=@zeros(6,6)
matrix staball=@zeros(2,11)
matrix deleteall=@zeros(216,7)
!delrow=0
matrix del
matrix paths
vector stabout
vector residuals
vector fitted
matrix infl
matrix stbexport
scalar numq
%stabgraphs=""
%inflgraphs=""
%delgraphs=""
for !which=1 to 2
  !m=!primary
  if !which=2 then
    !m=!trendmodel
  endif
  !run_model=!m
  !det=mgrid(!m,3)
  ord=@transpose(@subextract(mgrid,!m,6,!m,11))
  call bounds_boot(y,xbase,ord,!det,4,!B,79642,1,!m,boot)
  for !test=1 to 3
    !r=(!which-1)*3+!test
    iidout(!r,1)=!m
    iidout(!r,2)=!test
    for !j=1 to 4
      iidout(!r,!j+2)=boot(!test,!j)
    next
  next
  call design(y,xbase,ord,!det,4,1,w,v,nd)
  call olsfit(w,v,b,e,gi,ss,rk)
  call covariances(w,e,gi,ss,cv,hac,lev,cook)
  call influence_delete(w,v,nd,4,del)
  for !i=1 to @rows(del)
    !delrow=!delrow+1
    deleteall(!delrow,1)=!m
    for !j=1 to 6
      deleteall(!delrow,!j+1)=del(!i,!j)
    next
  next
  infl=@zeros(56,6)
  for !i=1 to 56
    infl(!i,1)=1968+!i
    infl(!i,2)=y(!i+4)
    infl(!i,3)=y(!i+4)-e(!i)
    infl(!i,4)=e(!i)
    infl(!i,5)=lev(!i)
    infl(!i,6)=cook(!i)
  next
  %iname="influence_"+@str(!which)
  table {%iname}
  %ifile="LCF_"+%iname+".csv"
  call export_table(infl,"year observed_lnLCF fitted_lnLCF residual leverage Cook_distance",%ifile,{%iname})
  smpl 1969 2024
  %resname="resid_"+@str(!which)
  %hatname="hat_"+@str(!which)
  %cookname="cook_"+@str(!which)
  series {%resname}
  series {%hatname}
  series {%cookname}
  mtos(e,{%resname})
  mtos(lev,{%hatname})
  mtos(cook,{%cookname})
  %g1="res_graph_"+@str(!which)
  %g2="hat_graph_"+@str(!which)
  %g3="cook_graph_"+@str(!which)
  freeze({%g1}) {%resname}.line
  freeze({%g2}) {%hatname}.line
  freeze({%g3}) {%cookname}.line
  %charttitle="Residuals: model "+@str(!m)
  {%g1}.addtext(t) %charttitle
  %charttitle="Leverage: model "+@str(!m)
  {%g2}.addtext(t) %charttitle
  %charttitle="Cook distance: model "+@str(!m)
  {%g3}.addtext(t) %charttitle
  %inflgraphs=%inflgraphs+" "+%g1+" "+%g2+" "+%g3
  for !width=1 to 5 step 4
    !dstart=1
    !dend=56
    if !width=5 then
      !dstart=57
      !dend=108
    endif
    vector dv=@subextract(del,!dstart,4,!dend,4)
    %dname="omit"+@str(!width)+"_"+@str(!which)
    !lastyear=2025-!width
    smpl 1969 {!lastyear}
    series {%dname}
    mtos(dv,{%dname})
    %dgname="delete_graph"+@str(!width)+"_"+@str(!which)
    freeze({%dgname}) {%dname}.line
    %charttitle="Interaction numerator; omitted block width "+@str(!width)+"; model "+@str(!m)
    {%dgname}.addtext(t) %charttitle
    %delgraphs=%delgraphs+" "+%dgname
  next
  call progress("Starting recursive stability",0,!BS,!run_model)
  call stability_boot(y,xbase,ord,!det,4,!BS,!m,stabout,paths)
  staball(!which,1)=!m
  for !j=1 to 10
    staball(!which,!j+1)=stabout(!j)
  next
  !firstyear=stabout(8)
  !sn=@rows(paths)
  stbexport=@zeros(!sn,7)
  for !i=1 to !sn
    stbexport(!i,1)=!firstyear+!i-1
    stbexport(!i,2)=paths(!i,1)
    stbexport(!i,3)=stabout(9)
    stbexport(!i,4)=-stabout(9)
    stbexport(!i,5)=paths(!i,2)
    stbexport(!i,6)=stabout(10)
    stbexport(!i,7)=-stabout(10)
  next
  %stname="stability_path_"+@str(!which)
  %stfile="LCF_"+%stname+".csv"
  table {%stname}
  call export_table(stbexport,"year CUSUM upper95 lower95 CUSUMSQ square_upper95 square_lower95",%stfile,{%stname})
  smpl {!firstyear} 2024
  for !stype=1 to 2
    !scol=2+3*(!stype=2)
    %prefix="st"+@str(!which)+"_"+@str(!stype)
    vector vpath=@columnextract(stbexport,!scol)
    vector vup=@columnextract(stbexport,!scol+1)
    vector vlo=@columnextract(stbexport,!scol+2)
    %vn=%prefix+"_path"
    %vu=%prefix+"_upper"
    %vl=%prefix+"_lower"
    series {%vn}
    series {%vu}
    series {%vl}
    mtos(vpath,{%vn})
    mtos(vup,{%vu})
    mtos(vlo,{%vl})
    %sgroup=%prefix+"_group"
    group {%sgroup} {%vn} {%vu} {%vl}
    %sg=%prefix+"_graph"
    freeze({%sg}) {%sgroup}.line
    %charttitle="Conditional bootstrap stability; model "+@str(!m)+"; type "+@str(!stype)
    {%sg}.addtext(t) %charttitle
    %stabgraphs=%stabgraphs+" "+%sg
  next
next
table iid_table
table stability_table
table deletion_table
call export_table(iidout,"model_id test_id statistic iid_p critical5 MCSE","LCF_iid_bootstrap.csv",iid_table)
call export_table(staball,"model_id supF split_year supF_p CUSUM CUSUM_p CUSUMSQ CUSUMSQ_p recursion_start CUSUM_critical CUSUMSQ_critical","LCF_stability.csv",stability_table)
call export_table(deleteall,"model_id width first_omitted last_omitted theta_inter theta_inter_p rho","LCF_deletion.csv",deletion_table)
results_spool.append iid_table stability_table deletion_table

' Exploratory marginal ECT-loading equations, HC3 covariance, conditional on gate.
table weak_table
if summary(!primary,24)=1 then
  ord=@transpose(@subextract(mgrid,!primary,6,!primary,11))
  call design(y,xbase,ord,1,4,1,w,v,nd)
  call olsfit(w,v,b,e,gi,ss,rk)
  vector(60) ect
  for !i=1 to 60
    ect(!i)=y(!i)+b(1)/b(2)
    for !j=1 to 5
      ect(!i)=ect(!i)+xbase(!i,!j)*b(!j+2)/b(2)
    next
  next
  matrix mw=@ones(56,8)
  vector(56) mv
  for !i=1 to 56
    !t=!i+4
    mw(!i,2)=ect(!t-1)
    for !j=1 to 3
      mw(!i,!j+2)=xbase(!t-1,!j)-xbase(!t-2,!j)
    next
    mw(!i,6)=y(!t-1)-y(!t-2)
    mw(!i,7)=xbase(!t,4)-xbase(!t-1,4)
    mw(!i,8)=(!t-1)/60
  next
  matrix weakout=@zeros(3,3)
  matrix meat
  matrix h3
  for !reg=1 to 3
    for !i=1 to 56
      mv(!i)=xbase(!i+4,!reg)-xbase(!i+3,!reg)
    next
    call olsfit(mw,mv,b,e,gi,ss,rk)
    meat=@zeros(8,8)
    for !i=1 to 56
      !hi=0
      for !j=1 to 8
        for !l=1 to 8
          !hi=!hi+mw(!i,!j)*gi(!j,!l)*mw(!i,!l)
        next
      next
      for !j=1 to 8
        for !l=1 to 8
          meat(!j,!l)=meat(!j,!l)+mw(!i,!j)*mw(!i,!l)*e(!i)^2/(1-!hi)^2
        next
      next
    next
    h3=gi*meat*gi
    weakout(!reg,1)=!reg
    weakout(!reg,2)=b(2)
    weakout(!reg,3)=2*(1-@ctdist(@abs(b(2)/@sqrt(h3(2,2))),48))
  next
  call export_table(weakout,"variable_id ECT_loading HC3_p","LCF_exploratory_exogeneity.csv",weak_table)
else
  weak_table(1,1)="Primary level gate failed; exploratory ECT-loading regressions not computed."
  weak_table.save(t=csv) LCF_exploratory_exogeneity.csv
endif
results_spool.append weak_table

call progress("Exporting figures and final report",0,1,!run_model)

' Annual-data figures: all graph objects are native and editable in EViews.
smpl @all
%seriesgraphs=""
for !j=1 to 6
  %vname=@word(%varnames,!j)
  %gname="series_graph_"+@str(!j)
  freeze({%gname}) {%vname}.line
  {%gname}.addtext(t) %vname
  %seriesgraphs=%seriesgraphs+" "+%gname
next
graph Fig01_series.merge {%seriesgraphs}
Fig01_series.align(2,0.6,0.6)
Fig01_series.save(t=png) LCF_Fig01_series.png
Fig01_series.save(t=pdf) LCF_Fig01_series.pdf
graph Fig04_influence.merge {%inflgraphs}
Fig04_influence.align(3,0.6,0.6)
Fig04_influence.save(t=png) LCF_Fig04_influence.png
Fig04_influence.save(t=pdf) LCF_Fig04_influence.pdf
graph Fig05_stability.merge {%stabgraphs}
Fig05_stability.align(2,0.6,0.6)
Fig05_stability.save(t=png) LCF_Fig05_stability.png
Fig05_stability.save(t=pdf) LCF_Fig05_stability.pdf
graph Fig06_deletion.merge {%delgraphs}
Fig06_deletion.align(2,0.6,0.6)
Fig06_deletion.save(t=png) LCF_Fig06_deletion.png
Fig06_deletion.save(t=pdf) LCF_Fig06_deletion.pdf
results_spool.append Fig01_series Fig04_influence Fig05_stability Fig06_deletion
wfsave LCF_NATIVE_RESULTS.wf1

' Preserve figure matrices before creating the unstructured comparison pages.
matrix fig_spec=@zeros(32,7)
for !i=1 to 32
  !critval=@qtdist(0.975,summary(!i,2)-summary(!i,3))
  fig_spec(!i,1)=summary(!i,8)
  fig_spec(!i,2)=summary(!i,8)-!critval*summary(!i,9)
  fig_spec(!i,3)=summary(!i,8)+!critval*summary(!i,9)
  fig_spec(!i,4)=summary(!i,15)
  fig_spec(!i,5)=summary(!i,16)
  fig_spec(!i,6)=summary(!i,17)
  fig_spec(!i,7)=0.05
next
matrix fig_coding=@zeros(8,6)
for !i=1 to 8
  !m=66+!i
  fig_coding(!i,1)=summary(!m,15)
  fig_coding(!i,2)=summary(!m,16)
  fig_coding(!i,3)=summary(!m,17)
  fig_coding(!i,4)=summary(!m,10)
  fig_coding(!i,5)=0.05
  fig_coding(!i,6)=summary(!m,8)
next
!gate=summary(!primary,24)
matrix fig_marginal=@zeros(3,6)
if !gate=1 then
  for !r=1 to !lrrow
    if lrall(!r,1)=!primary and lrall(!r,3)>=6 then
      !mr=lrall(!r,3)-5
      !mc=1+3*(lrall(!r,2)=2)
      fig_marginal(!mr,!mc)=lrall(!r,4)
      fig_marginal(!mr,!mc+1)=lrall(!r,6)
      fig_marginal(!mr,!mc+2)=lrall(!r,7)
    endif
  next
endif
' Metadata and final completion state are set only after all calculations succeed.
table run_information
run_information(1,1)="Software"
run_information(1,2)="EViews native only"
run_information(2,1)="Version number"
run_information(2,2)=@vernum
run_information(20,1)="Program build"
run_information(20,2)="2026-09-26-r6"
run_information(3,1)="B / B coding-stability / PSS simulations"
run_information(3,2)=@str(!B)+" / "+@str(!BS)+" / "+@str(!NPSS)
run_information(4,1)="Primary model ID"
run_information(4,2)=!primary
run_information(5,1)="Primary rule"
run_information(5,2)="BIC; L4; intercept; q0 allowed; original coding; common 1969-2024"
run_information(6,1)="Inference"
run_information(6,2)="Fixed-X; fixed selected model; separate F/t/Fx nulls; conditional association"
run_information(7,1)="No external estimation"
run_information(7,2)="No Python, R, COM calculation, shell computation or precomputed estimates"
run_information(8,1)="KPSS p-value bounds"
run_information(8,2)="p=0.01 denotes p<=0.01; p=0.10 denotes p>=0.10"
run_information(9,1)="Variable IDs"
run_information(9,2)="1 lnLCF; 2 lnGDPpc; 3 GFCF; 4 lnEnergy; 5 SI_z; 6 interaction"
run_information(10,1)="Long-run term IDs"
run_information(10,2)="1 GDP; 2 GFCF; 3 Energy; 4 SI; 5 interaction; 6-8 Energy marginal at SI=-1,0,1"
run_information(11,1)="Covariance IDs"
run_information(11,2)="1 OLS; 2 HAC2"
run_information(12,1)="Deterministic IDs"
run_information(12,2)="1 c; 2 ct; 3 c+post1979; 4 ct+post1979"
run_information(13,1)="Family IDs"
run_information(13,2)="1 main; 2 start1968; 3 qmin1; 4 ordinal coding"
run_information(14,1)="Coding IDs"
run_information(14,2)="0 original; 1 sqrt; 2 square; 3 log1p"
run_information(15,1)="Criterion IDs"
run_information(15,2)="1 BIC; 2 AICc"
run_information(16,1)="Bootstrap test IDs"
run_information(16,2)="1 F; 2 t; 3 Fx"
run_information(17,1)="Status"
run_information(17,2)="COMPUTATIONS COMPLETED; creating final comparison figures"
run_information.save(t=csv) LCF_run_information.csv

' Cross-page copies use explicit page-qualified names; no external data import.
pagecreate(page=specifications) u 32
copy annual\fig_spec fig_spec
copy annual\model_names model_names
series theta_inter
series lower95
series upper95
series pF
series pt
series pFx
series threshold
group specdata theta_inter lower95 upper95 pF pt pFx threshold
mtos(fig_spec,specdata)
group spec_inter theta_inter lower95 upper95
freeze(Fig02_interaction) spec_inter.line
Fig02_interaction.addtext(t) "Interaction level coefficient and 95% OLS limits; model IDs 1-32"
Fig02_interaction.save(t=png) LCF_Fig02_interaction.png
Fig02_interaction.save(t=pdf) LCF_Fig02_interaction.pdf
group spec_bounds pF pt pFx threshold
freeze(Fig03_bounds) spec_bounds.line
Fig03_bounds.addtext(t) "Three conditional bootstrap p-values; model IDs 1-32; threshold 0.05"
Fig03_bounds.save(t=png) LCF_Fig03_bounds.png
Fig03_bounds.save(t=pdf) LCF_Fig03_bounds.pdf
pagecreate(page=coding) u 8
copy annual\fig_coding fig_coding
series pF
series pt
series pFx
series inter_p
series threshold
series theta_inter
group codingdata pF pt pFx inter_p threshold theta_inter
mtos(fig_coding,codingdata)
group codingtests pF pt pFx inter_p threshold
freeze(Fig07_coding) codingtests.line
Fig07_coding.addtext(t) "Original, sqrt, square, log1p; each with c then ct; model IDs 67-74"
Fig07_coding.save(t=png) LCF_Fig07_coding.png
Fig07_coding.save(t=pdf) LCF_Fig07_coding.pdf
if !gate=1 then
  pagecreate(page=marginal) u 3
  copy annual\fig_marginal fig_marginal
  series ols_estimate
  series ols_lower
  series ols_upper
  series hac_estimate
  series hac_lower
  series hac_upper
  group marginaldata ols_estimate ols_lower ols_upper hac_estimate hac_lower hac_upper
  mtos(fig_marginal,marginaldata)
  freeze(Fig08_marginal) marginaldata.line
  Fig08_marginal.addtext(t) "Conditional long-run energy association; observations 1,2,3 denote SI_z=-1,0,+1"
  Fig08_marginal.save(t=png) LCF_Fig08_marginal.png
  Fig08_marginal.save(t=pdf) LCF_Fig08_marginal.pdf
endif
pageselect annual
copy specifications\Fig02_interaction annual\Fig02_interaction
copy specifications\Fig03_bounds annual\Fig03_bounds
copy coding\Fig07_coding annual\Fig07_coding
results_spool.append Fig02_interaction Fig03_bounds Fig07_coding
if !gate=1 then
  copy marginal\Fig08_marginal annual\Fig08_marginal
  results_spool.append Fig08_marginal
endif
run_information(17,2)="SUCCESS: all native computations and graph exports completed"
run_information(18,1)="Primary level gate passed"
run_information(18,2)=!gate
run_information(19,1)="Output directory"
run_information(19,2)=%outputdir
results_spool.append run_information
run_information.save(t=csv) LCF_run_information.csv
wfsave LCF_NATIVE_RESULTS.wf1
results_spool.save(t=rtf) LCF_results.rtf
show run_information
show results_spool

call progress("All computations and exports completed",1,1,!run_model)
run_status(1,1)="SUCCESS: all native calculations, workfile, tables and graphs saved."
run_status.save(t=csv) LCF_STATUS.csv
wfsave LCF_NATIVE_RESULTS.wf1
