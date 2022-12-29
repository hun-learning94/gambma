a=12; b=0.4; z=25; s=-10; nu = 1.1; theta = 0.14
u = seq(0, 1/nu, len=1e4)
samp = rtCCH(1e6, a, b, z, s, nu, theta)
hist(samp, nclass=100, probability = T,
     xlab= "u", ylab = "Density", col="#00c04b", border="white",
     main = paste0("tCCH(", a,", ", b,", ", z,", ", s,", ", nu,", ", theta, ")"))
lines(u, dtCCH(u, a, b, z, s, nu, theta), col="#008631", lwd=2, lty=2)
