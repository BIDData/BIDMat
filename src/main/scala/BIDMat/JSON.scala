
package BIDMat
import java.util.Map;
import com.cedarsoftware.util.io.JsonWriter;
import com.cedarsoftware.util.io.JsonReader;

object JSON {
  import java.util.Map;
  import java.util.HashMap; ;

  val jopts: Map[String, AnyRef] = new HashMap;
  jopts.put(JsonWriter.PRETTY_PRINT, java.lang.Boolean.TRUE);
  // jopts(JsonWriter.PRETTY_PRINT) = java.lang.Boolean.TRUE;

  def toJSON(obj: AnyRef): String = JsonWriter.objectToJson(obj);

  def toJSON(obj: AnyRef, pretty: Boolean): String = {
    if (pretty) {
      JsonWriter.objectToJson(obj, jopts);
    } else {
      JsonWriter.objectToJson(obj);
    }
  }

  def fromJSON(str: String): AnyRef = JsonReader.jsonToJava(str)
}