using EllipticFunctions
using SpecialFunctions
using Test

@testset "EllipticFunctions" begin
  @testset "Some values of the jtheta functions." begin
    @test isapprox(
      jtheta3(0, qfromtau(1im)), 
      pi^(1 / 4) / gamma(3 / 4)
    )
    @test isapprox(
      jtheta1(1 + 1im, qfromtau(1im)),
      1.1816128551455719 + 0.59589712760417439im
    )
    @test isapprox(
      jtheta1(20, qfromtau(0.01im)),
      0.03608696080206
    )
    @test isapprox(
      jtheta2(1 + 1im, qfromtau(1im)),
      0.74328632006610539 - 0.904159309718008im
    )
    @test isapprox(
      jtheta3(1 + 1im, qfromtau(1im)),
      0.86456184935441778 - 0.28488586703507289im
    )
    @test isapprox(
      jtheta4(1 + 1im, qfromtau(1im)),
      1.1351891564632007 + 0.28517396444192509im
    )
    @test real(jtheta1(1-1im, qfromtau(1.e-13*im))) == Inf
    @test imag(jtheta1(1-1im, qfromtau(1.e-13*im))) == -Inf
  end

  @testset "A value of jtheta1dash." begin
    @test isapprox(
      jtheta1dash(1 + 1im, qfromtau(1im)),
      0.81117649363854416 - 0.89452803853474627im
    )
  end

  @testset "Extended precision values of the jtheta functions." begin
    setprecision(512) do
      # Iex1 := N[EllipticTheta[3, 0, Exp[-Pi]], 160]
      # Iex1 = big"1.08643481121330801457531612151022345707020570724521888592079031598185673226710979596056162"
      @test isapprox(
        jtheta3(0, exp(-big(pi))), 
        pi^big"0.25" / gamma(big"0.75")
      )
      # Iex2 := N[EllipticTheta[1, 1 + I, Exp[-Pi], 160]
      Iex2 = (big"1.1816128551455718838220608128070906058460179459394417512927874046725942906428169906821482386001738910544957270663484565388200625773550622468047135600093911049900" +
        im*big"0.5958971276041743872161377011027469373635813856461054043464648167761758176537939322958996491529708881657121161701834193638614120732942080672217404420939827816559")
      @test isapprox(
        jtheta1(big"1" + 1im, big"1"*qfromtau(complex(0,big"1"))),
        Iex2)
      @test isapprox(jtheta1(big"1" + 1im, qfromtau(im)), Iex2)
      # Iex3 := N[EllipticTheta[1, 20, Exp[-\[Pi]/100]], 160]
      Iex3 = big"0.03608696080206363155776329635672833563464711788924314080325415263411276200554664115758304102825311403434396775483262414261489360017386937024719306449683254749470"
      @test isapprox(
        jtheta1(big"20", big"1"*qfromtau(im/big"100")),
        Iex3)
      # Iex4 := N[EllipticTheta[2, 1 + I, Exp[-Pi], 160]
      Iex4 = (big"0.7432863200661053824025802129786523834785362956592941951708805896822498191766850458788316332244978133281232831745818130600298973491832300423124975866346834329728" - 
        big"0.9041593097180079953651616403713868427230254526838112087993250064654302398913062980779292767687943979354760829763434496061036620743459953390408874339243536932842"*im)
      @test isapprox(
        jtheta2(big"1" + 1im, big"1"*qfromtau(complex(0,big"1"))),
        Iex4
      )
      # Iex5 := N[EllipticTheta[2, 1 + I, Exp[-Pi], 160]
      Iex5 = (big"0.8645618493544177805764400544500909422291535653000802605561669625295060578755985433578603618629183088323167111047150532663976491712264793925012387201665511406153" - 
        big"0.2848858670350728840140460807664851641159258489770164440718133953175273100784252022735901991395989448281649584146151824044349049550115567261773810228902173558426"*im)
      @test isapprox(
        jtheta3(big(1) + 1im, big"1"*qfromtau(complex(0,big"1"))),
        Iex5
      )
      # Iex6 := N[EllipticTheta[4, 1 + I, Exp[-Pi], 160]
      Iex6 = (big"1.1351891564632007223528568726888624612845443191876613951088991648371513534417278867794357360321619448328512180682186086252286983695970021076062865181006716567836" + 
        big"0.2851739644419250799934274072877530843224834859594948442834029203857299325079851140581878898547577578050570942789816372849270991723077875126594414684410369106887"*im)
      @test isapprox(
        jtheta4(big"1" + big"1"*im, exp(-big"1" * pi)),
        Iex6
      )
      @test isinf(real(jtheta1(big"1.0"-1im, qfromtau(eps(BigFloat)*im))))
      @test isinf(imag(jtheta1(big"1.0"-1im, qfromtau(eps(BigFloat)*im))))
      # Iex7 := N[EllipticThetaPrime[1, 1 + I, Exp[-Pi], 160]
      Iex7 = (big"0.8111764936385441957043074553366743391550522952772345100241804572577240584192926455803209936197529882979748406437548364116083492370624319996140934135203250303583" - 
        big"0.8945280385347462976675780747869855478742558275542709042707597767745675877012412975903384335407775407073407316973824753860860947118456849227212887278345194695573"*im)
      @test isapprox(
        jtheta1dash(big"1.0" + 1im, big"1"*qfromtau(big"1"*im)),
        Iex7
      )
    end
  end

  @testset "Periodicity-like properties of jtheta_ab." begin
    a   = 2 + 0.3im
    b   = 1 - 0.6im
    z   = 0.1 + 0.4im
    tau = 0.2 + 0.3im
    q   = qfromtau(tau)
    jab = jtheta_ab(a, b, z, q)
    @test isapprox(
      jtheta_ab(a, b, z + pi, q), 
      jab * exp(2im*a*pi)
    )
    @test isapprox(
      jtheta_ab(a, b, z + tau*pi, q), 
      jab * exp(-1im*(tau*pi + 2*z + 2*b*pi))
    )
  end

  @testset "Some values of etaDedekind." begin
    @test isapprox(etaDedekind(1im / 2), gamma(1 / 4) / 2^(7 / 8) / pi^(3 / 4))
    @test isapprox(etaDedekind(2im), gamma(1 / 4) / 2^(11 / 8) / pi^(3 / 4))
  end

  @testset "Lambda modular relation." begin
    x = 2.0
    @test isapprox(lambda(1im * sqrt(x)) + lambda(1im * sqrt(1 / x)), 1)
  end

  @testset "A value of kleinj." begin
    @test isapprox(
      kleinj(2im),
      66^3
    )
  end

  @testset "CarlsonRC and 'arc' functions." begin
    x = 1 + im
    y = 3
    @test isapprox(
      x * CarlsonRC(y*y, y*y-x*x),
      atanh(x/y)
    )
    @test isapprox(
      x * CarlsonRC(y*y-x*x, y*y),
      asin(x/y)
    )
    @test isapprox(
      x * CarlsonRC(y*y+x*x, y*y),
      asinh(x/y)
    )
    @test isapprox(
      sqrt(y*y-x*x) * CarlsonRC(x*x, y*y),
      acos(x/y)
    )
    @test isapprox(
      sqrt(x*x-y*y)*CarlsonRC(x*x, y*y),
      acosh(x/y)
    )
    y = -5 + 2im
    @test isapprox(
      x * CarlsonRC(y*y, y*y-x*x),
      -atanh(x/y)
    )
    @test isapprox(
      x * CarlsonRC(y*y-x*x, y*y),
      -asin(x/y)
    )
    @test isapprox(
      x * CarlsonRC(y*y+x*x, y*y),
      -asinh(x/y)
    )
    @test isapprox(
      sqrt(y*y-x*x) * CarlsonRC(x*x, y*y),
      acos(-x/y)
    )
    @test isapprox(
      sqrt(x*x-y*y) * CarlsonRC(x*x, y*y),
      acosh(-x/y)
    )
  end

  @testset "CarlsonRJ homogeneity." begin
    x = 1 + im
    y = -2 + 3im
    z = -3
    p = 4im
    kappa = 2
    @test isapprox(
      CarlsonRJ(x, y, z, p) / kappa / sqrt(kappa),
      CarlsonRJ(kappa*x, kappa*y, kappa*z, kappa*p)
    )
  end

  @testset "CarlsonRJ(x, y, y, p)." begin
    x = 1 + im
    y = -2 + 3im
    p = 4im
    @test isapprox(
      CarlsonRJ(x, y, y, p),
      3 * (CarlsonRC(x, y) - CarlsonRC(x, p)) / (p - y)
    )
  end

  @testset "A value of ellipticK." begin
    @test isapprox(
      ellipticK(0.5),
      8 * pi^(3/2) / gamma(-1/4)^2
    )
  end

  @testset "An extended precision value of ellipticK." begin
    @test isapprox(
      ellipticK(big"0.5"),
      8 * pi^(big"3"/2) / gamma(-big"1"/4)^2
    )
  end

  @testset "ellipticK and CarlsonRF." begin
    m = 2 - 3im
    @test isapprox(
      ellipticK(m),
      CarlsonRF(0, 1-m, 1)
    )
  end

  @testset "Complete ellipticE and CarlsonRG." begin
    m = 2 - 3im
    @test isapprox(
      ellipticE(m),
      2*CarlsonRG(0, 1-m, 1)
    )
  end

  @testset "Complete ellipticE and CarlsonRD." begin
    m = 2 - 3im
    @test isapprox(
      ellipticE(m),
      (1-m) * (CarlsonRD(0, 1-m, 1) + CarlsonRD(0, 1, 1-m)) / 3
    )
  end

  @testset "A value of complete ellipticE." begin
    @test isapprox(
      ellipticE(0.5),
      (2 * gamma(3/4)^4 + pi^2) / (4 * sqrt(pi) * gamma(3/4)^2)
    )
  end

  @testset "Symmetry of ellipticZ." begin
    phi = -5 + 3im
    m = -4 - 9im
    @test isapprox(
      ellipticZ(conj(phi), conj(m)),
      conj(ellipticZ(phi, m))
    )
  end

  @testset "Symmetry of ellipticZ." begin
    phi = -5 + 3im
    m = -4 - 9im
    @test isapprox(
      ellipticZ(conj(phi), conj(m)),
      conj(ellipticZ(phi, m))
    )
  end

  @testset "Relation ellipticZ/PI/E/K." begin
    phi = 7 - 6im
    m = -3
    @test isapprox(
      ellipticZ(phi, m),
      (1-m)*ellipticPI(phi, m, m) + m*sin(2*phi)/2/sqrt(1-m*sin(phi)^2) -
        ellipticE(m) / ellipticK(m) * ellipticPI(phi, 0, m)
    )
  end

  @testset "ellipticPI with n=1." begin
    phi = 1 + im
    m = 2 - im
    @test isapprox(
      ellipticPI(phi, 1, m),
      (sqrt(1-m*sin(phi)^2) * tan(phi) - ellipticE(phi,m)) / (1-m) +
        ellipticF(phi, m)
    )
  end

  @testset "Extended Precision ellipticF" begin
    setprecision(512) do
      phi = parse(BigFloat, "0.15")
      m = parse(BigFloat, "0.81")
      I = ellipticF(phi, m)
      IWolfram = parse(BigFloat, "0.15045731627390324557125914235016372659553548570625738854307684752660195708464381473254355631857978043355122247940117341312349152469759997869058851217919643931")
      @test abs(I - IWolfram) < 1e-154
    end
  end

  @testset "A value of agm." begin
    @test isapprox(
      agm(1, sqrt(2)),
      2 * pi^(3/2) * sqrt(2) / gamma(1/4)^2
    )
  end

  @testset "kleinjinv works." begin
    j = 0.1 + 0.2im
    @test isapprox(
      kleinj(kleinjinv(j)),
      j
    )
  end

  @testset "EisensteinE2 development." begin
    q = 0.005 + 0.005im
    @test isapprox(
      EisensteinE2(q),
      1 - 24 * (q/(1-q) + 2*q^2/(1-q^2) + 3*q^3/(1-q^3) + 4*q^4/(1-q^4) + 5*q^5/(1-q^5))
    )
  end

  @testset "Sum of the e_i is zero." begin
    omega1 = 1.4 - 1im
    omega2 = 1.6 + 0.5im
    omega = (omega1, omega2)
    e1 = wp(omega1; omega = omega)
    e2 = wp(omega2; omega = omega)
    e3 = wp(-omega1-omega2; omega = omega)
    @test isapprox(
      e1 + e2 + e3, 0; atol = 1e-10
    )
  end

  @testset "Differential equation wp." begin
    g2 = 1.4 - 1im
    g3 = 1.6 + 0.5im
    g = (g2, g3)
    z = 1 + 1im
    p = wp(z; g = g)
    pdash = wp(z; g = g, derivative = 1)
    @test isapprox(
      pdash^2,
      4 * p^3 - g2 * p - g3
    )
  end

  @testset "Equianharmonic wp." begin
    omega2 = gamma(1/3)^3 / 4 / pi
    z0     = omega2 * (1 + 1im/sqrt(3))
    @test isapprox(
      wp(z0; g = (0, 1)),
      0;
      atol = 10*eps(Float64), rtol = 0
    )
  end

  @testset "A value of wsigma." begin
    omega1 = gamma(1/4)^2 / 4 / sqrt(pi)
    omega2 = 1im * omega1
    omega = (omega1, omega2)
    @test isapprox(
      wsigma(omega1; omega = omega),
      exp(pi/8) * 2^(1/4)
    )
  end

  @testset "A value of wzeta given omega." begin
    omega1 = gamma(1/4)^2 / 4 / sqrt(pi)
    omega2 = 1im * omega1
    omega = (omega1, omega2)
    @test isapprox(
      wzeta(omega1; omega = omega),
      pi / 4 / omega1
    )
  end

  @testset "A value of wzeta given g." begin
    g = (5 + 3im, 5 + 3im)
    @test isapprox(
      wzeta(1 + 1im; g = g),
      0.802084165492408 - 0.381791358666872im
    )
  end

  @testset "Neville theta functions." begin
    z = 2.5
    m = 0.3
    @test isapprox(
      thetaC(z; m = m),
      -0.65900466676738154967
    )
    @test isapprox(
      thetaD(z; m = m),
      0.95182196661267561994
    )
    @test isapprox(
      thetaN(z; m = m),
      1.0526693354651613637
    )
    @test isapprox(
      thetaS(z; m = m),
      0.82086879524530400536
    )
  end

  @testset "Extended precision Neville theta functions." begin
    z = big"2.5"
    m = big"0.3"
    # Nex1 := N[NevilleThetaC[5/2, 3/10], 160]
    Nex1=big"-0.6590046667673815496485303765632134931323377681877931715936858355363291536494494386099152099197553390523718487512799789992274049195965925106214187963039501800221"
    @test isapprox(
      thetaC(z; m = m),
      Nex1
    )
    # Nex2 := N[NevilleThetaC[5/2, 3/10], 160]
    Nex2 = big"0.9518219666126756199502893873517518524498678861431149331777426589451505930972810721418127840410937816230582748951981263144179516821249111195948646242772204628267"
    @test isapprox(
      thetaD(z; m = m),
      Nex2
    )
    # Nex3 := N[NevilleThetaN[5/2, 3/10], 160]
    Nex3 = big"1.052669335465161363772100207227790438580960552534097530507833702184153288676597438817647515114228733873782994253724540324088812686517918938233601010706592302893"
    @test isapprox(
      thetaN(z; m = m),
      Nex3
    )
    # Nex4 := N[NevilleThetaS[5/2, 3/10], 160]
    Nex4 = big"0.8208687952453040055355508266515636852582693361520696044642407581282798088176868137965928997680175840294528934362468203055166725784495088747128187657877631004078"
    @test isapprox(
      thetaS(z; m = m),
      Nex4;
      atol = 1e-15
    )
  end

  @testset "jellip functions." begin
    u = 2 + 2im
    tau = 1im
    @test isapprox(
      jellip("cn", u; tau = tau)^2 + jellip("sn", u; tau = tau)^2,
      1
    )
  end

  @testset "am function." begin
    phi = 1 + 1im
    m = 2 - 1im
    u = ellipticF(phi, m)
    v = 2 - 2im
    psi = am(v, m)
    @test isapprox(
      am(u, m),
      phi
    )
    @test isapprox(
      ellipticF(psi, m),
      v
    )
  end

  @testset "Half-periods and elliptic invariants" begin
    x = gamma(1/4)^2 / 4 / sqrt(pi)
    g2, g3 = ellipticInvariants(x, im * x)
    @test isapprox(
      g2,
      1
    )
    @test isapprox(
      g3,
      0;
      atol = 1e-15
    )
    omega1, omega2 = halfPeriods(1, 0)
    g2, g3 = ellipticInvariants(omega1, omega2)
    @test isapprox(
      g2,
      1
    )
    @test isapprox(
      g3,
      0;
      atol = 1e-15
    )
  end
end
