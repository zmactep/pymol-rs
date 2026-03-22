const nav = document.querySelector("nav");
if (nav) {
  window.addEventListener("scroll", () => {
    nav.style.borderBottomColor =
      window.scrollY > 50 ? "var(--border-glow)" : "var(--border)";
  });
}
