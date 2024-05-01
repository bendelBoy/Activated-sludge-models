 {----------------------------------------------------------------------}
 {>>>>  DYNAMIC SIMULATION OF SINGLE CSTR ACTIVATED SLUDGE PROCESS  <<<<}
           {AUTHOR WILLI GUJER EAWAG DUEBENDORF SWITZERLAND}
 CONST IS  =3;                              {number of soluble compounds}
       IMAX=4;                                {total number of compounds}
       JMAX=2;                                      {number of processes}
 TYPE  tVECTOR  = ARRAY[1..IMAX]         OF REAL;
       tMATRIX  = ARRAY[1..IMAX]      OF tVECTOR;
       tNAME    = PACKED ARRAY[1..10]    OF CHAR;
 VAR   RHO    : ARRAY[1..JMAX]         OF REAL;           {process rates}
       NUE    : ARRAY[1..IMAX,1..JMAX] OF REAL;   {stoichiometric matrix}
       C0,C   : tVECTOR;                   {feed, effluent concentration}
       NAME   : ARRAY[1..IMAX]        OF tNAME;      {names of compounds}
       DH,DX,O2SAT,KLA,T : REAL;                      {system definition}
   {hydraulic-, solids-dilution rate, oxygen saturation, kla value, time}
 VAR   MUE,KS,KO,KN,B:REAL; {specific for problem to be solved: kinetics}
 {----------------------------------------------------------------------}
 {----------------------------------------------------------------------}
 PROCEDURE DEFINE_NAMES;              {SPECIFIC FOR PROBLEM TO BE SOLVED}
 BEGIN
   NAME[1]:='OXYGEN    ';        {>>>> OXYGEN MUST BE COMPOUND NO 1 <<<<}
   NAME[2]:='SUBSTRATE ';
   NAME[3]:='NURTRIENT ';
   NAME[4]:='BIOMASS   ';
 END;
 {----------------------------------------------------------------------}
 {----------------------------------------------------------------------}
 PROCEDURE DEFINE_STOICHIOMETRY;      {SPECIFIC FOR PROBLEM TO BE SOLVED}
 VAR   I,J   : INTEGER;       {running variable for compound and process}
       Y,INU : REAL;     {yield coefficient, nutrient content of biomass}
 BEGIN
   WRITE('YIELD                       = '); READLN(Y);
   WRITE('NUTRIENT CONTENT OF BIOMASS = '); READLN(INU);
   FOR I:=1 TO IMAX DO FOR J:=1 TO JMAX DO NUE[I,J]:=0;  {initialize nue}
   NUE[1,1]:=-(1-Y)/Y;
   NUE[2,1]:=-1/Y;
   NUE[3,1]:=-INU;
   NUE[4,1]:=1;
   NUE[1,2]:=-1;
   NUE[3,2]:=INU;
   NUE[4,2]:=-1;
 END; { OF DEFINE_STOICHIOMETRY }
 {----------------------------------------------------------------------}
 {----------------------------------------------------------------------}
 PROCEDURE DEFINE_KINETICS;           {SPECIFIC FOR PROBLEM TO BE SOLVED}
 BEGIN
   WRITE('MUEMAX  = '); READLN(MUE);
   WRITE('KS      = '); READLN(KS);
   WRITE('KO2     = '); READLN(KO);
   WRITE('KNUTR   = '); READLN(KN);
   WRITE('B DECAY = '); READLN(B);
 END; { OF DEFINE_STOICHIOMETRY }
 {----------------------------------------------------------------------}
 {----------------------------------------------------------------------}
 PROCEDURE GET_PROCESS_RATES;         {SPECIFIC FOR PROBLEM TO BE SOLVED}
 {>>>>   RHO MUST CONTAIN ALL COMPOUNDS WITH NEGATIVE NUE-VALUES    <<<<}
 BEGIN
   RHO[1]:=MUE*C[4]*(C[1]/(KO+C[1]))*(C[2]/(KS+C[2]))*C[3]/(KN+C[3]);
   RHO[2]:=B*C[4]*C[1]/(KO+C[1]);
 END; { OF GET_PROCESS_RATES }
 {----------------------------------------------------------------------}
 {----------------------------------------------------------------------}
 PROCEDURE DEFINE_REACTOR;                            {REMAINS UNCHANGED}
 BEGIN
   WRITE('HYDRAULIC DILUTION RATE     = '); READLN(DH);
   WRITE('SOLIDS DILUTION RATE        = '); READLN(DX);
   WRITE('KLA VALUE FOR AERATION      = '); READLN(KLA);
   WRITE('OXYGEN SATURATION VALUE     = '); READLN(O2SAT);
 END; { OF DEFINE_REACTOR }
 {----------------------------------------------------------------------}
 {----------------------------------------------------------------------}
 PROCEDURE DEFINE_FEED_CONCENTRATIONS;                {REMAINS UNCHANGED}
 CONST CMIN=1.0E-3;
     {minimum feed concentration for convergence of ITERATE_STEADY_STATE}
 VAR I         : INTEGER;
 BEGIN
   WRITELN('FEED CONCENTRATIONS');
   FOR I:=1 TO IMAX DO BEGIN
     WRITE(NAME[I]:10,'  C0(',I,')  = '); READLN(C0[I]);
     IF C0[I]<CMIN THEN C0[I]:=CMIN;
   END;
 END; { OF DEFINE_FEED_CONCENTRATIONS }
 {----------------------------------------------------------------------}
 {----------------------------------------------------------------------}
 PROCEDURE SOLVE(VAR D:tMATRIX;VAR F,R:tVECTOR);      {REMAINS UNCHANGED}
        {solves the linear equation system D*F=R with D[i,i]<>0}
 VAR I,K,J            : INTEGER;
     P,A              : REAL;
 BEGIN
   FOR I:=1 TO IMAX DO BEGIN
     P:=D[I,I];
     FOR K:=1 TO IMAX DO BEGIN
       A:=D[K,I]/P;
       IF K<>I THEN FOR J:=1 TO IMAX DO
         IF J<>I THEN D[K,J]:=D[K,J]-A*D[I,J];
     END;
     FOR K:=1 TO IMAX DO D[K,I]:=D[K,I]/P;
     FOR J:=1 TO IMAX DO D[I,J]:=D[I,J]/(-P);
     D[I,I]:=1/P;
   END;
   FOR I:=1 TO IMAX DO BEGIN
     R[I]:=0; FOR K:=1 TO IMAX DO R[I]:=R[I]+D[I,K]*F[K];
   END;
 END; { OF SOLVE }
 {----------------------------------------------------------------------}
 {----------------------------------------------------------------------}
 PROCEDURE ITERATE_STEADY_STATE;                      {REMAINS UNCHANGED}
 LABEL 1;
 CONST ITCRIT=1E-3;                    {criterion for iteration accuracy}
 VAR   I,I1,J     : INTEGER;
       SUM        : REAL;                            {iteration accuracy}
       K,KOLD     : ARRAY[1..JMAX] OF REAL;   {linearized rate constants}
       IL         : ARRAY[1..JMAX] OF INTEGER;       {limiting compounds}
       CONTROL    : tVECTOR;                    {input or control vector}
       SYSTEM     : tMATRIX;      {systems matrix: kinetics and effluent}
 BEGIN
   FOR J:=1 TO JMAX DO FOR I:=1 TO IMAX DO IF NUE[I,J]<0 THEN IL[J]:=I;
   FOR I:=1 TO IMAX DO CONTROL[I]:=-DH*C0[I];
   CONTROL[1]:=CONTROL[1]-KLA*O2SAT;           {oxygen input by aeration}
   FOR J:=1 TO JMAX DO K[J]:=1;    {initialize linearized rate constants}
   REPEAT  1:
     FOR I:=1 TO IMAX DO SYSTEM[1,I]:=0;      {initialize systems matrix}
     FOR I:=2 TO IMAX DO SYSTEM[I]:=SYSTEM[1];
     FOR I:=1 TO IMAX DO BEGIN                    {set up systems matrix}
       FOR J:=1 TO JMAX DO
         SYSTEM[I,IL[J]]:=SYSTEM[I,IL[J]]+K[J]*NUE[I,J];
       IF I<=IS THEN SYSTEM[I,I]:=SYSTEM[I,I]-DH      {transport soluble}
                ELSE SYSTEM[I,I]:=SYSTEM[I,I]-DX;      {transport solids}
     END;
     SYSTEM[1,1]:=SYSTEM[1,1]-KLA;                             {aeration}
     SOLVE(SYSTEM,CONTROL,C);
     FOR I:=1 TO IMAX DO IF C[I]<0 THEN    {check for limiting compounds}
       FOR J:=1 TO JMAX DO IF NUE[I,J]<0 THEN BEGIN
         IL[J]:=I; GOTO 1; END;      {begin with other limiting compound}
     GET_PROCESS_RATES;
     SUM:=0;  KOLD:=K;
     FOR J:=1 TO JMAX DO BEGIN
       K[J]:=RHO[J]/C[IL[J]];             {new linearized rate constants}
       SUM:=SUM+ABS(KOLD[J]/K[J]-1);
     END;
   UNTIL SUM<ITCRIT;
 END; { OF ITERATE_STEADY_STATE }
 {----------------------------------------------------------------------}
 {----------------------------------------------------------------------}
 PROCEDURE PRINT;                                     {REMAINS UNCHANGED}
 VAR I : INTEGER;
 BEGIN
   WRITELN('REACTOR CONCENTRATIONS ARE:');
   FOR I:=1 TO IMAX DO WRITELN(NAME[I]:10,'  C (',I,')  = ',C[I]:10:2);
 END; { OF PRINT }
 {----------------------------------------------------------------------}
 {----------------------------------------------------------------------}
 PROCEDURE INTEGRATION(TOLD,TNEW : REAL);             {REMAINS UNCHANGED}
 CONST MAXT1=0.02;{maximum change of single compound in time step at end}
       MAXT2=0.2;                                  {   "    at beginning}
       TDIF =0.1;                      {time to change to lower timestep}
 VAR I,J                : INTEGER;
     T,DELT,D,TO2,DTO2,DEL : REAL;
     R,CT               : tVECTOR;  {reaction rate, time derivative of C}
 BEGIN
   T:=TOLD;  DEL:=MAXT2;
   REPEAT                                              {integration loop}
     IF TNEW-T < TDIF THEN BEGIN DEL:=MAXT1; DELT:=TNEW-T END
                      ELSE DELT:=TNEW-T-0.99*TDIF;
     GET_PROCESS_RATES;
     FOR I:=2 TO IMAX DO BEGIN
       R[I]:=0; FOR J:=1 TO JMAX DO R[I]:=R[I]+RHO[J]*NUE[I,J];
       IF I>IS THEN D:=DX ELSE D:=DH;     {distinguish solid and soluble}
       CT[I]:=DH*C0[I]-D*C[I]+R[I];               {mass balance equation}
       IF C[I]<>0 THEN IF CT[I]<>0 THEN        {obtain allowed time step}
         IF ABS(DEL*C[I]/CT[I])<DELT THEN DELT:=ABS(DEL*C[I]/CT[I]);
     END;
     TO2:=T; T:=T+DELT;        {increment time, set time for oxygen loop}
     REPEAT                                 {integration loop for oxygen}
       DTO2:=T-TO2;
       GET_PROCESS_RATES;
       R[1]:=0; FOR J:=1 TO JMAX DO R[1]:=R[1]+NUE[1,J]*RHO[J];
       CT[1]:=DH*(C0[1]-C[1])+KLA*(O2SAT-C[1])+R[1];
       IF C[1]<>0 THEN IF CT[1]<>0 THEN        {obtain allowed time step}
         IF ABS(DEL*C[1]/CT[1])<DTO2 THEN DTO2:=ABS(DEL*C[1]/CT[1]);
       TO2:=TO2+DTO2; C[1]:=C[1]+CT[1]*DTO2;  {increment time and oxygen}
     UNTIL (TO2*1.00001>T);                             {rounding errors}
     FOR I:=2 TO IMAX DO C[I]:=C[I]+DELT*CT[I]; {increment concentration}
   UNTIL (T*1.000001>TNEW);
 END; { OF INTEGRATION }
 {----------------------------------------------------------------------}
 {----------------------------------------------------------------------}
 PROCEDURE DYNAMIC_SIMULATION;                        {REMAINS UNCHANGED}
 BEGIN
   REPEAT
     WRITE('TIME UNTIL NEXT OUTPUT      = '); READLN(T); WRITELN;
     IF T>0 THEN BEGIN                      {T=0 ends dynamic simulation}
       DEFINE_KINETICS; WRITELN;    {to adjust for change of temperature}
       DEFINE_REACTOR;  WRITELN;                     {to change flowrate}
       DEFINE_FEED_CONCENTRATIONS; WRITELN;  {change feed concentrations}
       INTEGRATION(0,T);
       PRINT; WRITELN;
     END;
   UNTIL T=0;                                       {to exit the program}
 END; { OF DYNAMIC_SIMULATION }
 {----------------------------------------------------------------------}
 {----------------------------------------------------------------------}
 BEGIN { OF MAIN PROGRAM }                            {REMAINS UNCHANGED}
   WRITELN('CHEMOSTAT SIMULATION PROGRAM'); WRITELN;
   DEFINE_NAMES;
   DEFINE_STOICHIOMETRY;                    WRITELN;
   DEFINE_KINETICS;                         WRITELN;
   DEFINE_REACTOR;                          WRITELN;
   DEFINE_FEED_CONCENTRATIONS;              WRITELN;
   ITERATE_STEADY_STATE;
   WRITELN('INITIAL CONDITIONS'); WRITELN('------------------');
   PRINT;                                   WRITELN;
   INTEGRATION(0,1/DX);        {integrate over one solids retention time}
   WRITELN('RESULT OF INTEGRATION'); WRITELN('--------------------');
   PRINT;                                   WRITELN;
   DYNAMIC_SIMULATION;
 END.  { OF MAIN PROGRAM }
 {----------------------------------------------------------------------}