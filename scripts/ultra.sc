
import com.pi4j.wiringpi.Gpio._

val trigPin = 4       
val echoPin = 5
val MAX_DISTANCE = 200          // define the maximum measured distance in cm 
val timeOut = MAX_DISTANCE * 60  // calculate timeout according to the maximum measured distance

def pulseIn(pin:Int, level:Int, timeout:Int):Int = 
{
  val m0 = micros();
  var t1 = 0L;
  var t2 = 0L;
  while (digitalRead(pin) != level && t1 < timeout) { 
    t1 = micros() - m0;
  }
  while (digitalRead(pin) == level && t2 < timeout) {
    t2 = micros() - m0;
  }
  if (t2 < timeout) { 
    t2.toInt
  } else { 
    0
  }
}

def getSonar():Double = {   //get the measurement result of ultrasonic module with unit: cm
  digitalWrite(trigPin,HIGH); //send 10us high level to trigPin 
  delayMicroseconds(10);
  digitalWrite(trigPin,LOW);
  val pingTime = pulseIn(echoPin,HIGH,timeOut);   //read plus time of echoPin
  val distance = pingTime * 340.0 / 2 / 10000.0; //calculate distance with sound speed 340m/s
  distance;
}

print("Program is starting ... \n");
wiringPiSetup();
pinMode(trigPin,OUTPUT);
pinMode(echoPin,INPUT);
while (true) {
  val distance = getSonar();
  print("The distance is : %.2f cm\n".format(distance));
  delay(1000);
}   
