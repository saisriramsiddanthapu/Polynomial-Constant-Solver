import java.io.*;
import java.math.BigInteger;
import java.util.*;

public class ConstantTerm {

  static class Pt { BigInteger x, y; Pt(BigInteger x, BigInteger y){ this.x=x; this.y=y; } }

  static String readAll() throws Exception {
    StringBuilder sb = new StringBuilder();
    try (BufferedReader br = new BufferedReader(new InputStreamReader(System.in))) {
      String line;
      while ((line = br.readLine()) != null) sb.append(line).append('\n');
    }
    return sb.toString();
  }

  static String unq(String s){
    s = s.trim();
    if (s.length()>=2 && s.charAt(0)=='"' && s.charAt(s.length()-1)=='"')
      return s.substring(1, s.length()-1);
    return s;
  }

  static class Parsed {
    int k;
    List<Pt> pts = new ArrayList<>();
  }

  static Parsed parse(String json) {
    String s = json;
    int keysIdx = s.indexOf("\"keys\"");
    int kStart = s.indexOf("\"k\"", keysIdx);
    int colon  = s.indexOf(':', kStart);
    int comma  = s.indexOf(',', colon);
    int close  = s.indexOf('}', colon);
    int endNum = (comma==-1?close:Math.min(comma, close));
    int k = Integer.parseInt(s.substring(colon+1, endNum).trim());

    List<Pt> pts = new ArrayList<>();
    int idx = 0;
    while (true){
      int q1 = s.indexOf('"', idx);
      if (q1==-1) break;
      int q2 = s.indexOf('"', q1+1);
      if (q2==-1) break;
      String key = s.substring(q1+1, q2);
      idx = q2+1;
      boolean numeric = key.chars().allMatch(Character::isDigit);
      if (!numeric) continue;
      int objStart = s.indexOf('{', idx);
      int objEnd   = -1;
      int depth=0;
      for (int i=objStart;i<s.length();i++){
        char c = s.charAt(i);
        if (c=='{') depth++;
        if (c=='}') { depth--; if (depth==0){ objEnd=i; break; } }
      }
      if (objStart==-1 || objEnd==-1) break;
      String obj = s.substring(objStart, objEnd+1);
      idx = objEnd+1;

      int bIdx = obj.indexOf("\"base\"");
      int bCol = obj.indexOf(':', bIdx);
      int bCom = obj.indexOf(',', bCol);
      String bStr = unq(obj.substring(bCol+1, bCom==-1?obj.length()-1:bCom).trim());

      int vIdx = obj.indexOf("\"value\"");
      int vCol = obj.indexOf(':', vIdx);
      int vEnd = obj.indexOf('}', vCol);
      String vStr = unq(obj.substring(vCol+1, vEnd).trim()).toLowerCase();

      BigInteger x = new BigInteger(key);
      int base = Integer.parseInt(bStr);
      BigInteger y = new BigInteger(vStr, base);
      pts.add(new Pt(x, y));
    }
    Parsed p = new Parsed(); p.k = k; p.pts = pts;
    return p;
  }

  static BigInteger lagrangeAtZero(List<Pt> pts, int k){
    LinkedHashMap<BigInteger, BigInteger> map = new LinkedHashMap<>();
    for (Pt p: pts) map.putIfAbsent(p.x, p.y);
    List<BigInteger> xs = new ArrayList<>(map.keySet()).subList(0, k);
    List<BigInteger> ys = new ArrayList<>();
    for (BigInteger x: xs) ys.add(map.get(x));

    BigInteger numSum = BigInteger.ZERO;
    BigInteger denLCM = BigInteger.ONE;

    List<BigInteger> termNums = new ArrayList<>();
    List<BigInteger> termDens = new ArrayList<>();

    for (int i=0;i<k;i++){
      BigInteger num = BigInteger.ONE;
      BigInteger den = BigInteger.ONE;
      for (int j=0;j<k;j++){
        if (i==j) continue;
        num = num.multiply(xs.get(j).negate());
        den = den.multiply(xs.get(i).subtract(xs.get(j)));
      }
      BigInteger tn = ys.get(i).multiply(num);
      BigInteger td = den;
      if (td.signum()<0){ td = td.negate(); tn = tn.negate(); }
      BigInteger g = tn.gcd(td);
      tn = tn.divide(g); td = td.divide(g);
      termNums.add(tn); termDens.add(td);
      denLCM = denLCM.multiply(td.divide(denLCM.gcd(td)));
    }
    for (int i=0;i<k;i++){
      BigInteger scale = denLCM.divide(termDens.get(i));
      numSum = numSum.add(termNums.get(i).multiply(scale));
    }
    BigInteger g = numSum.gcd(denLCM);
    numSum = numSum.divide(g); denLCM = denLCM.divide(g);
    if (!denLCM.equals(BigInteger.ONE))
      throw new RuntimeException("Non-integer constant");
    return numSum;
  }

  public static void main(String[] args) throws Exception {
    String json = readAll();
    Parsed p = parse(json);
    p.pts.sort(Comparator.comparing(a -> a.x));
    BigInteger c = lagrangeAtZero(p.pts, p.k);
    System.out.println(c.toString());
  }
}
